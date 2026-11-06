# =============================================================================
# Mouse scATAC-seq analysis
# =============================================================================
# Purpose:
#   Perform quality control, normalization, dimensionality reduction,
#   clustering, and gene-activity visualization for two mouse samples.
#
# Inputs:
#   - Aggregated 10x Genomics ATAC peak-barcode matrix and cell metadata.
#   - Fragment file and its accompanying tabix index.
#
# Project layout (all paths are relative to the project root):
#   data/atac/filtered_peak_bc_matrix.h5
#   data/atac/singlecell.csv
#   data/atac/fragments.tsv.gz
#   data/atac/fragments.tsv.gz.tbi
#   results/
#
# Outputs:
#   - results/atac_aggr.rds: annotated object before cell filtering.
#   - results/atac_subset.rds: filtered object with clustering and gene activity.
#
# Usage:
#   Run from the project root after configuring the paths and sample labels
#   below. Confirm the barcode-suffix order before assigning sample labels.
#   Quality-control thresholds and clustering settings are study-specific.
# =============================================================================

# 01. Load required packages ---------------------------------------------------

library(Signac)
library(Seurat)
library(hdf5r)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(clustree)

# 02. Configure input paths, output paths, and sample labels --------------------

# Edit these relative paths to match the local project layout.
input_dir <- file.path("data", "atac")
output_dir <- "results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Labels correspond to barcode suffixes in their first-occurrence order.
sample_labels <- c("Poly_Ai14+", "Poly_Ai14-")

# 03. Import the peak matrix and create the chromatin object --------------------

counts <- Read10X_h5(file.path(input_dir, "filtered_peak_bc_matrix.h5"))

metadata <- read.csv(
  file = file.path(input_dir, "singlecell.csv"),
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = file.path(input_dir, "fragments.tsv.gz"),
  min.cells = 1
)

atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  project = "ATAC",
  meta.data = metadata
)

# 04. Add gene annotations ----------------------------------------------------

# Retrieve mouse gene annotations and use UCSC-style chromosome names.
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
genome(annotations) <- "mm10"
Annotation(atac) <- annotations

# 05. Calculate quality-control metrics and assign sample labels ---------------

atac <- NucleosomeSignal(object = atac)
atac <- TSSEnrichment(atac)

# Extract the aggregation sample identifier from each cell-barcode suffix.
atac$sample_id <- gsub(".*-", "", colnames(atac))
current_samples <- unique(atac$sample_id)
print(current_samples)

# Review current_samples against the aggregation metadata before relabeling.
atac$sample_id <- factor(
  atac$sample_id,
  levels = current_samples,
  labels = sample_labels
)
table(atac$sample_id)

atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100
atac$blacklist_ratio <- atac$blacklist_region_fragments / atac$peak_region_fragments

VlnPlot(
  object = atac,
  features = c(
    "pct_reads_in_peaks",
    "peak_region_fragments",
    "TSS.enrichment",
    "blacklist_ratio",
    "nucleosome_signal"
  ),
  pt.size = 0.1,
  ncol = 5,
  group.by = "sample_id"
)

# Inspect empirical quantiles alongside the study-specific thresholds below.
nPercent_ReadsinPeaks <- quantile(
  atac$pct_reads_in_peaks,
  probs = c(0.01, 0.99)
)
nPercent_ReadsinPeaks

nPeak_Region_Fragments <- quantile(
  atac$peak_region_fragments,
  probs = c(0.05, 0.95)
)
nPeak_Region_Fragments

nTSS_enrichment <- quantile(atac$TSS.enrichment, probs = c(0.01, 0.99))
nTSS_enrichment

nNucleosome_Signal <- quantile(atac$nucleosome_signal, probs = c(0.05, 0.95))
nNucleosome_Signal

saveRDS(atac, file = file.path(output_dir, "atac_aggr.rds"), compress = TRUE)

# 06. Filter cells using the study-specific quality-control thresholds ---------

atac_subset <- subset(
  x = atac,
  subset =
    peak_region_fragments > 2600 &
    peak_region_fragments < 60000 &
    pct_reads_in_peaks > 72 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 1 &
    TSS.enrichment > 5
)

# 07. Normalize peak counts and perform dimensionality reduction ---------------

# Apply TF-IDF normalization and latent semantic indexing (LSI).
atac_subset <- RunTFIDF(atac_subset)
atac_subset <- FindTopFeatures(atac_subset, min.cutoff = "q0")
atac_subset <- RunSVD(atac_subset)
DepthCor(atac_subset)

# Use LSI components 2-30 for UMAP and nearest-neighbor graph construction.
atac_subset <- RunUMAP(atac_subset, reduction = "lsi", dims = 2:30)
atac_subset <- FindNeighbors(atac_subset, reduction = "lsi", dims = 2:30)

# 08. Evaluate clustering resolutions and inspect the selected clusters --------

atac_subset <- FindClusters(
  atac_subset,
  resolution = seq(0.1, 1.0, by = 0.1)
)
N_clustree <- clustree(atac_subset)
N_clustree

# Apply the resolution selected for the study dataset.
atac_subset <- FindClusters(atac_subset, resolution = 0.5, verbose = FALSE)
DimPlot(object = atac_subset, label = TRUE) + NoLegend()

p_state <- DimPlot(atac_subset, group.by = "sample_id", reduction = "umap")
p_cluster <- DimPlot(atac_subset, group.by = "seurat_clusters", label = TRUE)
p_state + theme(aspect.ratio = 1) | p_cluster + theme(aspect.ratio = 1)

# Save a checkpoint before adding the gene-activity assay.
saveRDS(atac_subset, file = file.path(output_dir, "atac_subset.rds"), compress = TRUE)

# 09. Calculate, normalize, and visualize gene activity ------------------------

# Gene activity estimates accessibility around genes; it is not measured RNA.
# Store these counts in an assay named "RNA" for downstream visualization.
gene.activities <- GeneActivity(atac_subset)
atac_subset[["RNA"]] <- CreateAssayObject(counts = gene.activities)
atac_subset <- NormalizeData(
  object = atac_subset,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = median(atac_subset$nCount_RNA)
)

DefaultAssay(atac_subset) <- "RNA"

FeaturePlot(
  object = atac_subset,
  features = c("Gzmk", "Eomes", "Tbx21", "Stat4"),
  order = TRUE,
  pt.size = 0.05,
  alpha = 0.8,
  max.cutoff = "q95",
  ncol = 2
)

# Update the checkpoint with the normalized gene-activity assay.
saveRDS(atac_subset, file = file.path(output_dir, "atac_subset.rds"), compress = TRUE)
