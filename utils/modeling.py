"""Repeated-split L1 logistic regression ("LASSO") models of CAR-T response.

Used by code/01_compute/01_ml_repeated_cv.ipynb and 02_ml_cohort3_holdout.ipynb.

Public API
----------
run_lasso_iterations(X, y, stratify, ...)
    Monte Carlo cross-validation: n_iter stratified 70/30 train/test splits.
    In each split, an L1 logistic regression is tuned by 5x5 repeated
    stratified CV on the training set and scored on the test set. Options:
    ``filter_75`` (expression filter), ``shuffle_labels`` (null control) and
    ``random_select30`` (random 30-gene control).
run_lasso_iterations_clinical(X, y, stratify, ...)
    Same splits, with a median-impute -> scale -> L1 logistic pipeline for
    clinical covariates.
run_lasso_iterations_LN_combined(X_clinical, y, stratify, X_LN, ...)
    Same splits. An out-of-fold LN-gene-set score plus clinical covariate(s)
    feed an L2 (default) or L1 logistic regression.
run_holdout_bootstrap(X, y, train_mask, test_mask, ...)
    A single model fitted on the train mask and evaluated on a disjoint test
    mask, with a bootstrap CI on the test AUC.
run_holdout_permutations(X, y, train_mask, test_mask, models, ...)
    Label-permutation null for run_holdout_bootstrap.
save_results(res, outdir)
    Writes the csv files read by the figure notebooks.

Reproducibility
---------------
Iteration ``i`` uses ``seed = i`` for the train/test split, the CV splitter,
the liblinear solver, the label shuffle (``seed`` and ``seed + 1``) and the
random gene choice. Every iteration is therefore independent of the others,
and ``iterations=range(5)`` gives exactly the first 5 of the 101 iterations.
The ``n_jobs`` setting has no effect on the results.
The splits are positional, so they also depend on
(a) the sample (row) order given by data_loading.load_modeling_data() and
(b) stratify labels with the same sorted order (add_stratify_group()).
"""
import os
from collections import namedtuple

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve
from sklearn.model_selection import (GridSearchCV, RepeatedStratifiedKFold,
                                     StratifiedKFold, train_test_split)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

# ---------------------------------------------------------------------------
# Shared constants and building blocks
# ---------------------------------------------------------------------------
C_GRID = np.logspace(-2, 2, 30)   # inverse regularisation strengths searched
CV_SPLITS = 5                     # inner CV: 5 folds ...
CV_REPEATS = 5                    # ... x 5 repeats
MAX_ITER = int(1e6)
POS_LABEL, NEG_LABEL = "R", "NR"

RepeatedCVResult = namedtuple("RepeatedCVResult", ["aucs", "classifier_scores", "selected_features"])
RepeatedCVResult.__doc__ = """Output of the run_lasso_iterations* functions.
Unpacks as an ``(auc_results, classifier_scores_df, selected_features_df)`` tuple."""

HoldoutResult = namedtuple("HoldoutResult", ["train_auc", "test_auc", "ci_low", "ci_high", "best_C",
                                             "scores_df", "features_df", "boot_aucs"])
HoldoutResult.__doc__ = """Output of run_holdout_bootstrap."""


def l1_logreg(seed, C=1.0):
    """L1 logistic regression (class-balanced liblinear) used in every model of the paper."""
    return LogisticRegression(penalty="l1", C=C, class_weight="balanced", max_iter=MAX_ITER,
                              solver="liblinear", random_state=seed)


def repeated_cv(seed):
    """5x5 repeated stratified K-fold splitter used to tune C."""
    return RepeatedStratifiedKFold(n_splits=CV_SPLITS, n_repeats=CV_REPEATS, random_state=seed)


def tune_C(estimator, X, y, seed, n_jobs, param="C"):
    """Grid-search ``param`` over C_GRID by ROC-AUC with repeated_cv(seed). Returns the fitted GridSearchCV."""
    gs = GridSearchCV(estimator, param_grid={param: C_GRID}, scoring="roc_auc",
                      cv=repeated_cv(seed), n_jobs=n_jobs)
    gs.fit(X, y)
    return gs


def filter_expressed(X_train, X_test, min_value=0.25, min_frac=0.25):
    """Keep the features with value > min_value in at least min_frac of the TRAIN samples (the 'filter_75' step)."""
    keep = (X_train > min_value).sum(axis=0) >= (min_frac * len(X_train))
    X_train = X_train.loc[:, keep]
    return X_train, X_test.loc[:, X_train.columns]


def youden_metrics(y_train, train_probs, y_test, test_probs):
    """Youden-optimal threshold chosen on TRAIN and applied unchanged to TEST.

    Returns (threshold, train_sens, train_spec, test_sens, test_spec).
    """
    fpr, tpr, thresholds = roc_curve(y_train, train_probs, pos_label=POS_LABEL)
    threshold = thresholds[np.argmax(tpr - fpr)]

    def sens_spec(y_true, probs):
        pred = np.where(probs >= threshold, POS_LABEL, NEG_LABEL)
        tn, fp, fn, tp = confusion_matrix(y_true, pred, labels=[NEG_LABEL, POS_LABEL]).ravel()
        sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
        spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
        return sens, spec

    return (threshold,) + sens_spec(y_train, train_probs) + sens_spec(y_test, test_probs)


# ---------------------------------------------------------------------------
# Repeated train/test splits (Monte Carlo CV)
# ---------------------------------------------------------------------------
def _align_stratify(stratify, X):
    if isinstance(stratify, pd.Series):
        return stratify.loc[X.index]
    return stratify


def _check_xy(X, y):
    if isinstance(y, pd.Series) and not y.index.equals(X.index):
        raise ValueError("y must be indexed like X (the split is positional)")


def _split(X, y, stratify_vec, test_size, seed, shuffle_labels):
    """Stratified train/test split for iteration ``seed``, with optional label shuffling."""
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, stratify=stratify_vec, random_state=seed)
    if shuffle_labels:
        # Permute label values but keep indices, so X/y alignment is preserved.
        y_train = pd.Series(y_train.sample(frac=1, random_state=seed).values, index=y_train.index)
        y_test = pd.Series(y_test.sample(frac=1, random_state=seed + 1).values, index=y_test.index)
    return X_train, X_test, y_train, y_test


def _repeated_splits(X, y, stratify, fit_iteration, n_iter, iterations, test_size,
                     shuffle_labels, add_youden, verbose):
    """Shared loop: split, call ``fit_iteration``, collect scores and coefficients.

    ``fit_iteration(X_train, X_test, y_train, seed)`` must return
    ``(train_probs, test_probs, feature_names, coefficients)``.
    """
    _check_xy(X, y)
    stratify_vec = _align_stratify(stratify, X)
    seeds = list(range(n_iter)) if iterations is None else list(iterations)

    aucs, scores, features = [], [], []
    for k, seed in enumerate(seeds):
        if verbose:
            print(f"Iteration {k + 1}/{len(seeds)} (seed {seed})...", end="\r")
        X_train, X_test, y_train, y_test = _split(X, y, stratify_vec, test_size, seed, shuffle_labels)
        train_probs, test_probs, feat_names, coefs = fit_iteration(X_train, X_test, y_train, seed)

        aucs.append(roc_auc_score(y_test, test_probs))

        df_train = {"sample_id": X_train.index, "iteration": seed, "train_or_test": "train",
                    "classifier_score": train_probs}
        df_test = {"sample_id": X_test.index, "iteration": seed, "train_or_test": "test",
                   "classifier_score": test_probs}
        if add_youden:
            thr, tr_sens, tr_spec, te_sens, te_spec = youden_metrics(y_train, train_probs, y_test, test_probs)
            df_train.update(threshold=thr, sensitivity=tr_sens, specificity=tr_spec)
            df_test.update(threshold=thr, sensitivity=te_sens, specificity=te_spec)
        scores.append(pd.concat([pd.DataFrame(df_train), pd.DataFrame(df_test)]))
        features.append(pd.DataFrame({"iteration": seed, "feature": feat_names, "coefficient": coefs}))

    if verbose:
        print(f"\nDone: {len(seeds)} iterations, median test AUC = {np.median(aucs):.3f}")
    return RepeatedCVResult(aucs, pd.concat(scores, ignore_index=True),
                            pd.concat(features, ignore_index=True))


def _nonzero(names, coefs):
    coefs = np.asarray(coefs).flatten()
    mask = coefs != 0
    return np.asarray(names)[mask], coefs[mask]


def run_lasso_iterations(X, y, stratify, n_iter=101, iterations=None, test_size=0.3,
                         filter_75=False, shuffle_labels=False, random_select30=False,
                         n_jobs=30, verbose=True):
    """Repeated stratified train/test splits with a CV-tuned L1 logistic regression.

    Parameters
    ----------
    X : DataFrame, samples x features (normalised expression).
    y : Series of "R"/"NR", indexed like X.
    stratify : Series (aligned to X by index) or array used to stratify the split.
        In the paper this is ``response + "_" + cohort``.
    n_iter : number of splits; iteration i uses seed i.
    iterations : optional iterable of seeds to run instead of range(n_iter).
    test_size : test fraction (0.3).
    filter_75 : keep features with value > 0.25 in >= 25% of train samples
        (computed on the training split only).
    shuffle_labels : permute y_train and y_test after the split (null control).
    random_select30 : after filtering, keep 30 random features (control),
        drawn with RandomState(seed).
    n_jobs : GridSearchCV workers (has no effect on the results).

    Returns
    -------
    RepeatedCVResult(aucs, classifier_scores, selected_features). The
    classifier_scores columns are sample_id, iteration, train_or_test,
    classifier_score, threshold, sensitivity, specificity. The threshold is the
    Youden cut-off from the training split.
    """
    def fit_iteration(X_train, X_test, y_train, seed):
        if filter_75:
            X_train, X_test = filter_expressed(X_train, X_test)
        else:
            # Keep this copy/column selection: the memory layout of the frames
            # changes the last bits of predict_proba (BLAS), so it is kept for
            # identical scores.
            X_train = X_train.copy()
            X_test = X_test.loc[:, X_train.columns]
        if random_select30:
            # Same draws as np.random.seed(seed); np.random.choice(...).
            genes = np.random.RandomState(seed).choice(X_train.columns, size=30, replace=False)
            X_train, X_test = X_train.loc[:, genes], X_test.loc[:, genes]
        best = tune_C(l1_logreg(seed), X_train, y_train, seed, n_jobs).best_estimator_
        names, coefs = _nonzero(X_train.columns, best.coef_)
        return best.predict_proba(X_train)[:, 1], best.predict_proba(X_test)[:, 1], names, coefs

    return _repeated_splits(X, y, stratify, fit_iteration, n_iter, iterations, test_size,
                            shuffle_labels, add_youden=True, verbose=verbose)


def clinical_pipeline(seed, penalty="l1"):
    """Median-impute -> StandardScaler -> logistic regression (L1/liblinear or L2/lbfgs).

    The last step is named 'lasso' (L1) or 'logreg' (L2).
    """
    if penalty == "l1":
        step = ("lasso", l1_logreg(seed))
    elif penalty == "l2":
        step = ("logreg", LogisticRegression(penalty="l2", class_weight="balanced", max_iter=MAX_ITER,
                                             solver="lbfgs", random_state=seed))
    else:
        raise ValueError("penalty must be 'l1' or 'l2'")
    return Pipeline([("imputer", SimpleImputer(strategy="median")),
                     ("scaler", StandardScaler()),
                     step])


def run_lasso_iterations_clinical(X, y, stratify, n_iter=101, iterations=None, test_size=0.3,
                                  shuffle_labels=False, n_jobs=20, verbose=True):
    """run_lasso_iterations for clinical covariates, using clinical_pipeline(penalty='l1').

    Imputation and scaling are fitted inside each CV fold, so there is no
    leakage. classifier_scores has no Youden columns.
    """
    def fit_iteration(X_train, X_test, y_train, seed):
        best = tune_C(clinical_pipeline(seed, "l1"), X_train, y_train, seed, n_jobs,
                      param="lasso__C").best_estimator_
        names, coefs = _nonzero(X.columns, best.named_steps["lasso"].coef_)
        return best.predict_proba(X_train)[:, 1], best.predict_proba(X_test)[:, 1], names, coefs

    return _repeated_splits(X, y, stratify, fit_iteration, n_iter, iterations, test_size,
                            shuffle_labels, add_youden=False, verbose=verbose)


def oof_gene_score(X_train, y_train, X_test, seed, n_jobs, n_oofs=5):
    """LN-type gene-set score for a combined model.

    1. Tune C for an L1 model on the full training split (5x5 CV).
    2. Out-of-fold train scores: StratifiedKFold(n_oofs, shuffle, seed) with that C.
    3. Test scores: the model refitted on the whole training split with that C.
    Returns (oof_train_scores, test_scores), both P(R).
    """
    best_C = tune_C(l1_logreg(seed), X_train, y_train, seed, n_jobs).best_params_["C"]

    X_arr, y_arr = X_train.to_numpy(), y_train.to_numpy()
    oof = np.full(len(y_train), np.nan)
    for fold_tr, fold_val in StratifiedKFold(n_splits=n_oofs, shuffle=True,
                                             random_state=seed).split(X_arr, y_arr):
        m = l1_logreg(seed, C=best_C).fit(X_arr[fold_tr], y_arr[fold_tr])
        oof[fold_val] = m.predict_proba(X_arr[fold_val])[:, 1]   # classes_ = [NR, R]

    full = l1_logreg(seed, C=best_C).fit(X_train, y_train)
    return oof, full.predict_proba(X_test)[:, 1]


def run_lasso_iterations_LN_combined(X_clinical, y, stratify, X_LN, n_iter=101, iterations=None,
                                     test_size=0.3, n_oofs=5, final_penalty="l2", n_jobs=20,
                                     verbose=True):
    """Combined model of [out-of-fold LN score | clinical covariate(s)].

    The split is made on X_clinical with the same seeds as the other analyses.
    The LN score comes from oof_gene_score(). The combined matrix then goes
    through clinical_pipeline(final_penalty), tuned by 5x5 CV.
    * final_penalty='l2' (results/ml/forced_clin_ridge, SFig10B): ridge. Every
      coefficient is saved, so the clinical variable is always kept.
    * final_penalty='l1': LASSO. Only the nonzero coefficients are saved.
    """
    step = "lasso" if final_penalty == "l1" else "logreg"

    def fit_iteration(X_train_clin, X_test_clin, y_train, seed):
        oof, test_ln = oof_gene_score(X_LN.loc[X_train_clin.index], y_train,
                                      X_LN.loc[X_test_clin.index], seed, n_jobs, n_oofs)
        X_train = pd.concat([pd.Series(oof, index=X_train_clin.index, name="LN_score").to_frame(),
                             X_train_clin], axis=1)
        X_test = pd.concat([pd.Series(test_ln, index=X_test_clin.index, name="LN_score").to_frame(),
                            X_test_clin], axis=1)
        best = tune_C(clinical_pipeline(seed, final_penalty), X_train, y_train, seed, n_jobs,
                      param=step + "__C").best_estimator_
        coefs = best.named_steps[step].coef_.flatten()
        names = list(X_train.columns)
        if final_penalty == "l1":
            names, coefs = _nonzero(names, coefs)
        return best.predict_proba(X_train)[:, 1], best.predict_proba(X_test)[:, 1], names, coefs

    return _repeated_splits(X_clinical, y, stratify, fit_iteration, n_iter, iterations, test_size,
                            shuffle_labels=False, add_youden=False, verbose=verbose)


# ---------------------------------------------------------------------------
# Single train/test holdout with bootstrap CI
# ---------------------------------------------------------------------------
def _as_mask(mask, X):
    if isinstance(mask, pd.Series):
        mask = mask.loc[X.index]
    mask = np.asarray(mask, dtype=bool)
    if mask.shape != (len(X),):
        raise ValueError("mask must have one entry per row of X")
    return mask


def bootstrap_auc(y_true, probs, n_bootstrap, random_state=42):
    """Bootstrap the AUC of fixed predictions. Draws with a single class are skipped."""
    rng = np.random.RandomState(random_state)
    y_arr, p_arr = np.asarray(y_true), np.asarray(probs)
    boot = []
    for _ in range(n_bootstrap):
        idx = rng.randint(0, len(y_arr), len(y_arr))
        if len(np.unique(y_arr[idx])) < 2:
            continue
        boot.append(roc_auc_score(y_arr[idx], p_arr[idx]))
    return boot


def run_holdout_bootstrap(X, y, train_mask, test_mask, gene_list=None, filter_75=False,
                          n_bootstrap=1000, random_state=42, n_jobs=20,
                          train_label="train", test_label="test", verbose=True):
    """Fit one L1 logistic model on ``train_mask`` and evaluate it on ``test_mask``.

    C is tuned by 5x5 repeated stratified CV on the training samples only
    (seed ``random_state``). The model is fitted once. Only the test
    evaluation is bootstrapped (``n_bootstrap`` resamples, RandomState
    ``random_state``) to give a 95% percentile CI on the test AUC.

    The masks are boolean Series (aligned to X by index) or arrays, and they
    must be disjoint. In the paper: train = Cohort 1 + 2, test = Cohort 3.
    """
    train_mask, test_mask = _as_mask(train_mask, X), _as_mask(test_mask, X)
    assert not np.any(train_mask & test_mask), "train and test masks overlap"
    assert train_mask.any() and test_mask.any(), "empty train or test set"

    X_train, y_train = X.loc[train_mask], y.loc[train_mask]
    X_test, y_test = X.loc[test_mask], y.loc[test_mask]
    if verbose:
        print(f"Train ({train_label}): {len(y_train)} samples  |  Test ({test_label}): {len(y_test)} samples")

    if gene_list is not None:
        cols = X.columns.intersection(gene_list)
        X_train, X_test = X_train[cols], X_test[cols]
        if verbose:
            print(f"Gene list: {len(gene_list)} requested, {len(cols)} found in data")
    if filter_75:
        X_train, X_test = filter_expressed(X_train, X_test)
        if verbose:
            print(f"After filter_75: {X_train.shape[1]} features retained")

    gs = tune_C(l1_logreg(random_state), X_train, y_train, random_state, n_jobs)
    best, best_C = gs.best_estimator_, gs.best_params_["C"]

    train_probs = best.predict_proba(X_train)[:, 1]
    test_probs = best.predict_proba(X_test)[:, 1]
    train_auc = roc_auc_score(y_train, train_probs)
    test_auc = roc_auc_score(y_test, test_probs)

    boot = bootstrap_auc(y_test, test_probs, n_bootstrap, random_state)
    ci_low, ci_high = np.percentile(boot, 2.5), np.percentile(boot, 97.5)
    if verbose:
        print(f"Best C: {best_C:.4f} | Train AUC ({train_label}): {train_auc:.3f} | "
              f"Test AUC ({test_label}): {test_auc:.3f}, 95% CI {ci_low:.3f}-{ci_high:.3f} "
              f"(n={len(boot)} valid resamples)")

    scores_df = pd.concat([
        pd.DataFrame({"sample_id": X_train.index, "split": "train", "classifier_score": train_probs}),
        pd.DataFrame({"sample_id": X_test.index, "split": "test", "classifier_score": test_probs}),
    ], ignore_index=True)
    names, coefs = _nonzero(X_train.columns, best.coef_)
    features_df = pd.DataFrame({"feature": names, "coefficient": coefs})
    return HoldoutResult(train_auc, test_auc, ci_low, ci_high, best_C, scores_df, features_df, boot)


def run_holdout_permutations(X, y, train_mask, test_mask, models, n_perm=100, n_bootstrap=200,
                             random_state=42, outdir=None, n_jobs=20, verbose=False):
    """Label-permutation null for run_holdout_bootstrap.

    For each of ``n_perm`` permutations, all labels are permuted with one
    RandomState(random_state) stream, and then
    every model in ``models`` is fitted in order. ``models`` maps a name to
    run_holdout_bootstrap keyword arguments, e.g.
    ``{"LN": {"gene_list": LN_ENSEMBL}, "LASSO": {"filter_75": True}}``.
    If ``outdir`` is given, each fit is saved to ``outdir/perm_<name>/iter_<i>/``.
    Returns a DataFrame with the columns perm, <name>_test_auc, ...
    """
    rng = np.random.RandomState(random_state)
    rows = []
    for i in range(n_perm):
        y_perm = pd.Series(rng.permutation(y.values), index=y.index, name=y.name)
        row = {"perm": i}
        for name, kwargs in models.items():
            res = run_holdout_bootstrap(X, y_perm, train_mask, test_mask, n_bootstrap=n_bootstrap,
                                        n_jobs=n_jobs, verbose=verbose, **kwargs)
            row[f"{name}_test_auc"] = res.test_auc
            if outdir is not None:
                save_results(res, os.path.join(outdir, f"perm_{name}", f"iter_{i}"))
        rows.append(row)
        if (i + 1) % 10 == 0:
            print(f"Completed {i + 1}/{n_perm} permutations")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def save_results(res, outdir):
    """Write a result to ``outdir`` (the file names and columns read by the figure notebooks).

    RepeatedCVResult -> classifier_scores.csv, selected_features.csv
    HoldoutResult    -> classifier_scores.csv, selected_features.csv,
                        bootstrap_aucs.csv, auc_summary.csv
    """
    os.makedirs(outdir, exist_ok=True)
    if isinstance(res, HoldoutResult):
        res.scores_df.to_csv(os.path.join(outdir, "classifier_scores.csv"), index=False)
        res.features_df.to_csv(os.path.join(outdir, "selected_features.csv"), index=False)
        pd.DataFrame({"bootstrap_auc": res.boot_aucs}).to_csv(
            os.path.join(outdir, "bootstrap_aucs.csv"), index=False)
        pd.DataFrame([{"train_auc": res.train_auc, "test_auc": res.test_auc, "ci_low": res.ci_low,
                       "ci_high": res.ci_high, "best_C": res.best_C}]).to_csv(
            os.path.join(outdir, "auc_summary.csv"), index=False)
    elif isinstance(res, RepeatedCVResult):
        res.classifier_scores.to_csv(os.path.join(outdir, "classifier_scores.csv"), index=False)
        res.selected_features.to_csv(os.path.join(outdir, "selected_features.csv"), index=False)
    else:
        raise TypeError(f"unsupported result type: {type(res)}")
