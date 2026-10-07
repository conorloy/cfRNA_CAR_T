# cfRNA_CAR_T

Code for Loy *et al.*, "Cell-free RNA Signatures Derived from the Tumor Microenvironment Predict Outcomes of CAR-T
Therapy in Large B Cell Lymphoma", *Cell Reports Medicine* (2026).

## Layout

```
config.yaml     paths to data/, results/ and plots/
environment/    conda environments
utils/          shared R and Python functions
code/
  00_scrna_prep/   download and subset the public scRNA-seq atlases
  01_compute/      ML models, differential abundance, clinical models -> results/
  02_figures/      figure notebooks -> plots/
data/           cfRNA counts and sample metadata
results/        outputs of code/01_compute (ML results included)
```

## Setup

```bash
conda env create -f environment/environment_R.yml         # R notebooks, kernel "ir"
conda env create -f environment/environment_modeling.yml  # ML notebooks, kernel "scvi"
conda env create -f environment/environment_scanpy.yml    # scRNA notebooks, kernel "tabsap"
```

Register each environment as a Jupyter kernel with the name in quotes, e.g.
`Rscript -e 'IRkernel::installspec(name = "ir")'` and `python -m ipykernel install --user --name scvi`.

## Data

The cfRNA count matrices and sample metadata are in `data/`. Raw sequencing data are available from GEO
([GSE325315](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325315)). The scRNA-seq atlases used for the dot plots
(Tabula Sapiens; Li *et al.* 2025; Strati *et al.* 2023) are downloaded from CELLxGENE by
`code/00_scrna_prep/00_download_public_data.ipynb`.

## Running

Run notebooks from their own folder, e.g.

```bash
cd code/01_compute
jupyter nbconvert --to notebook --execute --inplace 03_differential_abundance.ipynb
```

1. `00_scrna_prep/`: `00_download_public_data`, then `01_subset_dotplot_atlases` (only needed for `Fig3_dotplots`).
2. `01_compute/`: `01_ml_repeated_cv` and `02_ml_cohort3_holdout` (results included; ~4 h to rerun),
   `03_differential_abundance`, `04_clinical_adjustment_mice`.
3. `02_figures/`: any order.

| Figure | Notebook |
|---|---|
| Fig. 1, S1 | `Fig1_cohort_overview` |
| Fig. 2, S2–S6 | `Fig2_LN_module` |
| Fig. 3A–B, S3B–C, S7, S8, S10–S12 | `Fig3_ML_performance` |
| Fig. 3C, S13 | `Fig3_dotplots` |
| Fig. 4, S14–S19 | `Fig4_differential_abundance` |
| S9, Tables S3–S4 | `SFig9_clinical_adjustment`, `01_compute/04_clinical_adjustment_mice` |

## Citation

Loy C, Agun G, Maurer K, Bosch Vilaseca A, Potapova D, Jacobson CA, Ritz J, De Vlaminck I. Cell-free RNA Signatures
Derived from the Tumor Microenvironment Predict Outcomes of CAR-T Therapy in Large B Cell Lymphoma.
*Cell Reports Medicine* (2026).

## License

MIT
