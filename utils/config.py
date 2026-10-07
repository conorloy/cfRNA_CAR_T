"""Path configuration shared by all Python notebooks.

Reads <repo>/config.yaml, then the optional, untracked <repo>/config.local.yaml, whose keys
override those of config.yaml. Exposes in_path(), res_path(), res_in() and plot_path().
"""
from pathlib import Path

import yaml


def find_repo_root(start=None):
    d = Path(start or Path.cwd()).resolve()
    for p in [d, *d.parents]:
        if (p / "config.yaml").exists():
            return p
    raise FileNotFoundError(f"config.yaml not found above {d}")


REPO_ROOT = find_repo_root()


def _resolve(p):
    p = Path(p).expanduser()
    return (p if p.is_absolute() else REPO_ROOT / p).resolve()


with open(REPO_ROOT / "config.yaml") as fh:
    _cfg = yaml.safe_load(fh)
if (REPO_ROOT / "config.local.yaml").exists():
    with open(REPO_ROOT / "config.local.yaml") as fh:
        _cfg.update(yaml.safe_load(fh) or {})

INPUT_DIR = _resolve(_cfg["input_dir"])
RESULTS_DIR = _resolve(_cfg["results_dir"])
PLOTS_DIR = _resolve(_cfg["plots_dir"])

if not INPUT_DIR.exists():
    print(f"WARNING: input_dir does not exist: {INPUT_DIR} (edit config.yaml or config.local.yaml)")


def in_path(*parts):
    """File in the input data directory; errors early if it is missing."""
    p = INPUT_DIR.joinpath(*parts)
    if not p.exists():
        raise FileNotFoundError(f"Input file not found: {p}")
    return p


def res_path(*parts):
    """File in the results directory, for writing; parent folders are created on demand."""
    p = RESULTS_DIR.joinpath(*parts)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def res_in(*parts):
    """Existing file or folder in the results directory, for reading; errors if it is missing."""
    p = RESULTS_DIR.joinpath(*parts)
    if not p.exists():
        raise FileNotFoundError(f"Result not found: {p} (run the code/01_compute or "
                                "code/00_scrna_prep notebook that writes it)")
    return p


def plot_path(*parts):
    """File in the plots directory."""
    p = PLOTS_DIR.joinpath(*parts)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p
