"""Human CD8+ T-cell atlas integration with scVI.

Purpose
-------
Combine ten standardized CD8+ T-cell cohorts and export an integrated AnnData
object and trained model for downstream analysis.

Inputs
------
One quality-controlled CD8+ AnnData file per GEO dataset, with human gene symbols
as var_names, raw counts in layers['counts'] (or X), and harmonized cell metadata.
Cohort-specific cell selection and gene-identifier harmonization precede this
script. Input sample_id values must already distinguish samples across cohorts.

Project layout (relative to the project root)
--------------------------------------------
data/human_cd8/GSE137029.h5ad ... GSE297584.h5ad
results/human_cd8/

Usage
-----
Configure the relative paths below, then run:
    python Python_script_for_human_CD8_integration.py
See README_human_CD8_integration.md for the input contract and reference setup.

Outputs
-------
Human_CD8_scVI_integrated.h5ad: counts, metadata and obsm['X_scVI'].
scvi_model/: trained scVI model, loadable with the exported AnnData.
integration_summary.json: dimensions, parameters and software versions.
"""

# 01. Load required packages --------------------------------------------------

import json
import platform
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import skmisc  # Required by Scanpy's seurat_v3 feature-selection method.
import torch
from scipy import sparse


# 02. Configure input paths and study parameters ------------------------------

INPUT_DIR = Path("data") / "human_cd8"
OUTPUT_DIR = Path("results") / "human_cd8"

# Preserve the cohort order used in the original integration.
DATASET_IDS = (
    "GSE137029", "GSE145926", "GSE206299", "GSE283471", "GSE289496",
    "GSE131907", "GSE132465", "GSE176078", "GSE193816", "GSE297584",
)

# GSE146264 was excluded because its available matrix was not comparable to
# the UMI-scale count matrices used in the primary atlas.
PROJECT_GENES = (
    "GZMK", "EOMES", "CD8A", "CD8B", "GZMB", "PRF1", "NKG7",
    "TCF7", "CCR7", "IL7R", "TOX", "PDCD1", "CXCR6",
)
REQUIRED_OBS = (
    "dataset_id", "donor_id", "sample_id", "library_id", "disease", "tissue",
)
SEED = 42
MIN_SHARED_GENES = 7000
N_TOP_GENES = 4000
MIN_CELLS_PER_SAMPLE = 10
MODEL_PARAMETERS = {
    "n_hidden": 256,
    "n_latent": 30,
    "n_layers": 2,
    "gene_likelihood": "nb",
}
TRAINING_PARAMETERS = {
    "max_epochs": 200,
    "batch_size": 512,
    "early_stopping": True,
    "early_stopping_patience": 20,
    "check_val_every_n_epoch": 1,
}


# 03. Read count-backed cohorts and check identifiers -------------------------

def check_counts(counts, label, seed=SEED, sample_size=1_000_000):
    """Check count values without transforming or rounding the input matrix."""
    values = counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
    if values.size == 0:
        raise ValueError(f"{label}: empty count matrix.")
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError(f"{label}: counts must be finite and nonnegative.")
    rng = np.random.default_rng(seed)
    sampled = values[rng.choice(values.size, min(values.size, sample_size), False)]
    if np.mean(np.isclose(sampled, np.round(sampled))) < 0.99:
        raise ValueError(f"{label}: counts fail the integer-like input check.")


def read_cohort(dataset_id, path):
    """Read one standardized CD8+ cohort and retain counts and metadata."""
    source = sc.read_h5ad(path)
    missing = [column for column in REQUIRED_OBS if column not in source.obs]
    if missing:
        raise ValueError(f"{dataset_id}: missing metadata columns {missing}.")
    if source.obs[list(REQUIRED_OBS)].isna().any().any():
        raise ValueError(f"{dataset_id}: required metadata contains missing values.")
    if set(source.obs["dataset_id"].astype(str)) != {dataset_id}:
        raise ValueError(f"{dataset_id}: unexpected dataset_id values.")
    if not source.obs_names.is_unique:
        raise ValueError(f"{dataset_id}: cell identifiers are not unique.")
    if not source.var_names.is_unique:
        source.var_names_make_unique()

    # X is an allowed fallback only when it already contains unnormalized counts.
    counts = source.layers["counts"] if "counts" in source.layers else source.X
    check_counts(counts, dataset_id)
    cohort = ad.AnnData(
        X=counts.copy(), obs=source.obs.copy(), var=source.var.copy()
    )
    cohort.layers["counts"] = cohort.X.copy()

    if "orig_cell_id" not in cohort.obs:
        cohort.obs["orig_cell_id"] = cohort.obs_names.astype(str)
    prefix = f"{dataset_id}__"
    cohort.obs_names = [
        str(cell) if str(cell).startswith(prefix) else f"{prefix}{cell}"
        for cell in cohort.obs_names
    ]
    if not cohort.obs_names.is_unique:
        raise ValueError(f"{dataset_id}: duplicate cell identifiers after prefixing.")
    return cohort


# 04. Merge cohorts on their shared gene space --------------------------------

def assemble_atlas(input_dir=INPUT_DIR):
    """Concatenate cohorts in the original cell and gene order."""
    cohorts = [
        read_cohort(dataset_id, input_dir / f"{dataset_id}.h5ad")
        for dataset_id in DATASET_IDS
    ]
    shared = set.intersection(*(set(cohort.var_names) for cohort in cohorts))
    shared_ordered = [gene for gene in cohorts[0].var_names if gene in shared]
    if len(shared_ordered) < MIN_SHARED_GENES:
        raise ValueError(f"Only {len(shared_ordered)} genes are shared across cohorts.")
    missing = sorted(set(PROJECT_GENES) - shared)
    if missing:
        raise ValueError(f"Project genes missing from the shared gene space: {missing}")

    atlas = ad.concat(
        [cohort[:, shared_ordered].copy() for cohort in cohorts],
        join="inner", merge="same", index_unique=None,
    )
    if not atlas.obs_names.is_unique or not atlas.var_names.is_unique:
        raise ValueError("Merged cell and gene identifiers must be unique.")
    sample_cohorts = atlas.obs.groupby("sample_id", observed=True)["dataset_id"].nunique()
    if (sample_cohorts > 1).any():
        raise ValueError("sample_id values must distinguish samples across cohorts.")
    check_counts(atlas.layers["counts"], "Merged atlas", seed=123, sample_size=2_000_000)
    return atlas


# 05. Select integration features before removing small samples ---------------

def select_integration_features(atlas):
    """Select dataset-aware HVGs and retain the prespecified project genes."""
    atlas.X = atlas.layers["counts"].copy()
    sc.pp.highly_variable_genes(
        atlas,
        n_top_genes=N_TOP_GENES,
        flavor="seurat_v3",
        batch_key="dataset_id",
        layer="counts",
        subset=False,
    )
    atlas.var["force_retained_project_gene"] = atlas.var_names.isin(PROJECT_GENES)
    atlas.var.loc[list(PROJECT_GENES), "highly_variable"] = True
    # The union can exceed N_TOP_GENES if a project gene was not already an HVG.
    selected = atlas.var["highly_variable"].astype(bool)
    return atlas[:, selected].copy()


def filter_small_samples(atlas):
    """Remove samples with fewer than ten retained CD8+ cells after HVG selection."""
    sample_counts = atlas.obs["sample_id"].astype(str).value_counts()
    retained = sample_counts[sample_counts >= MIN_CELLS_PER_SAMPLE].index
    filtered = atlas[atlas.obs["sample_id"].isin(retained)].copy()
    if filtered.n_obs == 0:
        raise ValueError("No cells remain after the sample-size filter.")
    filtered.uns["sample_filter"] = {
        "min_cells_per_sample": MIN_CELLS_PER_SAMPLE,
        "n_samples_removed": int((sample_counts < MIN_CELLS_PER_SAMPLE).sum()),
        "n_cells_removed": int(atlas.n_obs - filtered.n_obs),
    }
    for column in ("dataset_id", "sample_id", "donor_id", "disease", "tissue"):
        filtered.obs[column] = pd.Categorical(filtered.obs[column].astype(str))
    return filtered


# 06. Train scVI and obtain the integrated representation ----------------------

def fit_scvi(atlas):
    """Train on raw counts using sample_id as the batch variable."""
    scvi.settings.seed = SEED
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(SEED)
    torch.set_float32_matmul_precision("high")

    # Disease and tissue remain biological metadata, without covariate correction.
    scvi.model.SCVI.setup_anndata(atlas, layer="counts", batch_key="sample_id")
    model = scvi.model.SCVI(atlas, **MODEL_PARAMETERS)
    model.train(**TRAINING_PARAMETERS)
    atlas.obsm["X_scVI"] = model.get_latent_representation()
    return model


# 07. Export the integrated object, model, and concise provenance --------------

def main():
    """Run the study integration through a downstream-ready AnnData object."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    atlas = assemble_atlas()
    n_shared_genes = atlas.n_vars
    atlas = select_integration_features(atlas)
    atlas = filter_small_samples(atlas)
    model = fit_scvi(atlas)

    software = {
        "python": platform.python_version(),
        **{
            package: version(package)
            for package in (
                "anndata", "numpy", "pandas", "scanpy", "scipy", "scvi-tools",
                "scikit-misc", "torch",
            )
        },
    }
    summary = {
        "n_cells": int(atlas.n_obs),
        "n_genes": int(atlas.n_vars),
        "n_shared_genes_before_hvg": int(n_shared_genes),
        "n_datasets": int(atlas.obs["dataset_id"].nunique()),
        "n_samples": int(atlas.obs["sample_id"].nunique()),
        "n_dataset_specific_donors": int(
            atlas.obs[["dataset_id", "donor_id"]].drop_duplicates().shape[0]
        ),
        "input_datasets": list(DATASET_IDS),
        "excluded_dataset": "GSE146264",
        "feature_selection": {
            "method": "seurat_v3", "n_top_genes": N_TOP_GENES,
            "batch_key": "dataset_id", "project_genes": list(PROJECT_GENES),
            "order": "shared genes, HVGs plus project genes, sample-size filter",
        },
        "sample_filter": atlas.uns["sample_filter"],
        "scvi": {
            "layer": "counts", "batch_key": "sample_id", "seed": SEED,
            "model_parameters": MODEL_PARAMETERS,
            "training_parameters": TRAINING_PARAMETERS,
        },
        "software": software,
    }
    atlas.uns["integration"] = summary
    model.save(OUTPUT_DIR / "scvi_model", overwrite=True)
    atlas.write_h5ad(OUTPUT_DIR / "Human_CD8_scVI_integrated.h5ad", compression="gzip")
    with (OUTPUT_DIR / "integration_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    print(f"Saved {atlas.n_obs:,} cells x {atlas.n_vars:,} genes to {OUTPUT_DIR}.")


if __name__ == "__main__":
    main()
