# Multi-omics analyses of GZMK-expressing CD8+ T cells

This code collection contains the mouse single-cell RNA-seq, paired TCR V(D)J,
single-cell ATAC-seq and cross-modality analyses used to characterize
GZMK-expressing CD8+ T cells. A separate Python workflow integrates public human
CD8+ T-cell datasets with scVI and exports an object for downstream analysis.

## Software and system requirements

Mouse single-cell RNA-seq data were processed upstream using Cell Ranger
(version 9.0.1) and analyzed in R using Seurat (version 5.4.0). The RNA workflow
performs quality control, clustering, cell-type annotation, marker analysis and
paired TCR clonotype assignment from Cell Ranger V(D)J annotations.

Mouse single-cell ATAC-seq data were processed upstream using Cell Ranger ATAC
(version 2.2.0) and analyzed using Signac (version 1.16.0) and Seurat. The ATAC
workflow includes quality control, TF-IDF normalization, latent semantic indexing,
clustering and gene-activity analysis. Cross-modality analysis uses Seurat label
transfer, Signac coverage and differential-accessibility analysis, and chromVAR
motif deviations with JASPAR2020 motifs and the mouse mm10 reference genome.

The human CD8+ T-cell integration reference environment used Python 3.11.15,
scanpy 1.11.5, anndata 0.12.19, scvi-tools 1.4.2, numpy 2.4.6, pandas 2.3.3,
scipy 1.17.1, scikit-misc 0.5.2 and PyTorch 2.13.0. The recorded environment ran
on Linux with an NVIDIA RTX 5090 GPU and a CUDA 13.0 PyTorch build. Use a PyTorch
build compatible with the target hardware.

The R software versions above are recorded in the study methods; Python versions
are from the recorded integration environment. The exact R version and versions
of additional R packages are not specified here. All required package imports
are listed at the beginning of each script. The R dependencies are:

- **CRAN:** Seurat, Signac, hdf5r, dplyr, readr, tibble, ggplot2, patchwork,
  scales, clustree, pheatmap and rlang.
- **Bioconductor:** EnsDb.Mmusculus.v79, GenomeInfoDb, GenomicRanges, JASPAR2020,
  TFBSTools, BSgenome.Mmusculus.UCSC.mm10, motifmatchr and chromVAR.

Install the packages required by the selected workflow before execution. The RNA
script requires a Seurat version supporting `JoinLayers()` and the `layer`
argument. Record `sessionInfo()` for each R analysis to preserve the actual
execution environment. The Python workflow records its package versions in
`integration_summary.json`.

## Scripts

| Script | Purpose |
| --- | --- |
| [R_script_for_scRNA_analyses.R](R_script_for_scRNA_analyses.R) | Mouse scRNA-seq analysis and paired TCR clonotype assignment. |
| [R_script_for_scATAC_analyses.R](R_script_for_scATAC_analyses.R) | Mouse scATAC-seq preprocessing, clustering and gene-activity analysis. |
| [R_script_for_scRNA_and_scATAC_intergrated_analyses.R](R_script_for_scRNA_and_scATAC_intergrated_analyses.R) | RNA-to-ATAC label transfer, coverage visualization, motif analysis and comparison of differential accessibility with RNA expression. |
| [Python_script_for_human_CD8_integration.py](Python_script_for_human_CD8_integration.py) | Integration of ten human CD8+ T-cell cohorts with scVI, ending with export of the integrated AnnData object and trained model. |

## Input data

Use this directory as the project root, or edit the relative paths in Section 02
of each script. Prepared input data are not included in this code collection.

```text
data/
  gene_expression/Sample1/ ... Sample4/   # 10x feature-barcode matrices
  vdj/Sample1/ ... Sample4/               # filtered_contig_annotations.csv
  atac/
    filtered_peak_bc_matrix.h5
    singlecell.csv
    fragments.tsv.gz
    fragments.tsv.gz.tbi
  objects/
    Poly_RNA.rds
    Poly_ATAC.rds
  human_cd8/
    <GEO_accession>.h5ad
results/
```

For mouse RNA analysis, the default sample mapping is `Sample1` = 2 DPI Ai14
positive, `Sample2` = 2 DPI Ai14 negative, `Sample3` = 14 DPI Ai14 positive and
`Sample4` = 14 DPI Ai14 negative. DPI denotes days after the last OVA inhalation.
Confirm the mapping against the prepared inputs. The ATAC script uses two sample
labels, `Poly_Ai14+` and `Poly_Ai14-`, assigned from aggregation barcode suffixes.

The cross-modality script requires preprocessed Seurat objects. The RNA object
must contain an `RNA` assay and `sample`/`seurat_clusters` metadata. The ATAC
object must contain a `peaks` chromatin assay with mm10 gene annotations, LSI and
UMAP reductions, and `sample_id`/`seurat_clusters` metadata. Its fragment file
and tabix index must be available. Check the objects and configure their paths:
the RNA and ATAC preprocessing scripts do not automatically create the filenames
`Poly_RNA.rds` and `Poly_ATAC.rds`.

The human workflow expects one standardized, quality-controlled CD8+ AnnData file
for each of GSE137029, GSE145926, GSE206299, GSE283471, GSE289496, GSE131907,
GSE132465, GSE176078, GSE193816 and GSE297584. Each input must contain harmonized
human gene symbols, raw counts in `layers['counts']` (or raw counts in `X` if
that layer is absent), and the `obs` columns `dataset_id`, `donor_id`, `sample_id`,
`library_id`, `disease` and `tissue`. Sample identifiers must distinguish samples
across cohorts. Cell selection and gene-identifier harmonization precede this
script.

## Usage

Run the R scripts section by section in R or RStudio, in their numbered order.
Review the configured sample labels, cluster selections and thresholds before
execution. The cross-modality workflow uses the prepared mouse RNA and ATAC
objects; human atlas integration is an independent workflow.

Run the human integration script from the project root:

```bash
python Python_script_for_human_CD8_integration.py
```

This workflow selects dataset-aware highly variable genes, retains prespecified
project genes, removes samples with fewer than ten retained cells, and trains
scVI using raw counts with `sample_id` as the batch variable. It exports a
30-dimensional latent representation for downstream analysis. UMAP, clustering
and downstream figure generation are outside this Python script's scope.

## Outputs

- **Mouse RNA:** `results/Poly_top10_cluster_markers.csv`,
  `results/Poly_CD8_markers.csv` and
  `results/ab_T_merge_selected_genes.rds`. The RDS is a clustered checkpoint
  saved before cell-type annotation and TCR metadata assignment. The final
  TCR-annotated object remains in the R session and requires an explicit
  `saveRDS()` call if it is to be retained.
- **Mouse ATAC:** `results/atac_aggr.rds` and `results/atac_subset.rds`.
  The latter contains the filtered, clustered object and normalized gene activity.
- **Mouse cross-modality analysis:** updated Seurat objects, differential-analysis
  tables and plots remain in the R session; the script does not export them.
- **Human integration:** `results/human_cd8/Human_CD8_scVI_integrated.h5ad`,
  `results/human_cd8/scvi_model/` and
  `results/human_cd8/integration_summary.json`. The AnnData object retains raw
  counts for the selected integration genes and stores the integrated
  representation in `obsm['X_scVI']`; its `X` matrix contains raw counts.
