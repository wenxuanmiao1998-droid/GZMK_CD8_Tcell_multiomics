# =============================================================================
# Mouse scRNA-seq and scATAC-seq integration and comparative analysis
# =============================================================================
# Purpose:
#   Transfer RNA-derived cell-type labels to ATAC cells using gene activity,
#   visualize chromatin accessibility, calculate motif activity with chromVAR,
#   and compare differential peak accessibility with RNA expression.
#
# Inputs:
#   - Preprocessed Seurat objects containing mouse scRNA-seq and scATAC-seq data.
#   - An ATAC fragment file and its accompanying tabix index.
#   - Mouse mm10 gene annotations in the ATAC object and the mm10 genome package.
#
# Project layout (all paths are relative to the project root):
#   data/objects/Poly_ATAC.rds
#   data/objects/Poly_RNA.rds
#   data/atac/fragments.tsv.gz
#   data/atac/fragments.tsv.gz.tbi
#
# Outputs:
#   Objects, marker tables, peak annotations, and plots remain in the R session.
#   This script does not write analysis results to disk.
#
# Usage:
#   Run from the project root after configuring the input paths below.
#   The ATAC object must contain existing LSI and UMAP reductions, sample_id and
#   seurat_clusters metadata, and a "peaks" chromatin assay with gene annotation.
#   The RNA object must contain an RNA assay and sample/seurat_clusters metadata.
#   Sample selections, cluster IDs, cell-type labels, and thresholds are
#   study-specific. Review these settings before reusing this workflow.
# =============================================================================

# 01. Load required packages ---------------------------------------------------

library(Seurat)
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)
library(chromVAR)
library(GenomeInfoDb)
library(GenomicRanges)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(scales)
library(pheatmap)
library(rlang)

# 02. Configure input paths ---------------------------------------------------

# Edit these relative paths to match the local project layout.
object_dir <- file.path("data", "objects")
fragment_file <- file.path("data", "atac", "fragments.tsv.gz")

# 03. Import preprocessed objects and reconnect ATAC fragments -----------------

ATAC <- readRDS(file.path(object_dir, "Poly_ATAC.rds"))
RNA <- readRDS(file.path(object_dir, "Poly_RNA.rds"))

# Use the chromatin assay for fragment connections and gene activity.
DefaultAssay(ATAC) <- "peaks"

# Replace the stored fragment connection with the configured local file.
# Cell names must match the barcodes represented by the fragment object.
fragments_new <- CreateFragmentObject(
  path = fragment_file,
  cells = colnames(ATAC)
)
Fragments(ATAC) <- NULL
Fragments(ATAC) <- fragments_new

# 04. Inspect existing ATAC clusters and sample proportions --------------------

# 04.1. Plot existing UMAP embeddings with labels at cluster medians.
My_Panel <- c(
  "#a2d2e7", "#67abcd", "#3581b7", "#b3e19b", "#6fb3a8",
  "#50aa4b", "#cdb6da", "#704ba3", "#9a7fdb", "#dba9a8"
)
names(My_Panel) <- levels(Idents(ATAC))
emb <- Embeddings(ATAC, "umap") %>%
  as.data.frame()
emb$cluster <- factor(Idents(ATAC))
centers <- emb %>%
  group_by(cluster) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )

PDIM <- DimPlot(ATAC, label = TRUE, pt.size = 0.01, cols = My_Panel)
PDIM +
  geom_point(
    data = centers,
    aes(umap_1, umap_2, colour = cluster),
    shape = 21,
    fill = alpha("white", 0.5),
    size = 6,
    stroke = 0.3,
    show.legend = FALSE
  ) +
  geom_text(
    data = centers,
    aes(umap_1, umap_2, label = cluster),
    color = "black",
    size = 4
  )

PDIM_SEP <- DimPlot(
  ATAC,
  pt.size = 0.01,
  cols = My_Panel,
  split.by = "sample_id"
) &
  theme(legend.position = "none")
PDIM_SEP +
  geom_point(
    data = centers,
    aes(umap_1, umap_2, colour = cluster),
    shape = 21,
    fill = alpha("white", 0.5),
    size = 6,
    stroke = 0.3,
    show.legend = FALSE
  ) +
  geom_text(
    data = centers,
    aes(umap_1, umap_2, label = cluster),
    color = "black",
    size = 4
  )

# 04.2. Calculate the proportion of cells in each cluster within each sample.
meta <- ATAC@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(cluster = Idents(ATAC))
df_prop <- meta %>%
  group_by(sample_id, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(freq = n / sum(n))

ggplot(df_prop, aes(x = sample_id, y = freq, fill = cluster)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Fraction of cells") +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = My_Panel) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 15)
  )

# 05. Select the RNA reference and identify cross-modality anchors -------------

# Restrict the RNA reference to the two 14-day sample groups.
Idents(RNA) <- "sample"
RNA <- subset(RNA, idents = c("14dpi_Ai14pos", "14dpi_Ai14neg"))
RNA <- FindVariableFeatures(RNA, nfeatures = 3000)

# Gene activity estimates accessibility around genes; it is not measured RNA.
gene.activities <- GeneActivity(ATAC, features = VariableFeatures(RNA))
ATAC[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(ATAC) <- "ACTIVITY"
ATAC <- NormalizeData(ATAC)
ATAC <- ScaleData(ATAC, features = rownames(ATAC))
transfer.anchors <- FindTransferAnchors(
  reference = RNA,
  query = ATAC,
  features = VariableFeatures(RNA),
  reference.assay = "RNA",
  query.assay = "ACTIVITY",
  reduction = "cca"
)

# 06. Annotate RNA clusters and transfer cell-type labels to ATAC cells --------

# Labels correspond to the current RNA cluster-level order.
celltype <- c(
  "central memory T",
  "exhausted T(1)",
  "memory T(1)",
  "GZMK+ inflammatory T",
  "memory T(2)",
  "GZMK+ stem-like T",
  "GZMK+ memory T",
  "exhausted T(2)"
)
Idents(RNA) <- "seurat_clusters"
names(celltype) <- levels(RNA)
RNA <- RenameIdents(RNA, celltype)
RNA$celltype <- Idents(RNA)
ref_label <- RNA$celltype

# Weight transferred labels using ATAC LSI components 2-30.
pred <- TransferData(
  anchorset = transfer.anchors,
  refdata = ref_label,
  weight.reduction = ATAC[["lsi"]],
  dims = 2:30
)
ATAC <- AddMetaData(ATAC, metadata = pred)
cells_use <- WhichCells(
  ATAC,
  expression = predicted.id %in% c(
    "exhausted T(1)", "memory T(1)", "memory T(2)", "GZMK+ memory T"
  )
)
My_Panel2 <- c(
  "exhausted T(1)" = "#67abcd",
  "memory T(1)" = "#3581b7",
  "memory T(2)" = "#b3e19b",
  "GZMK+ memory T" = "#ffc17f"
)
PPRDICT <- DimPlot(
  ATAC,
  cells = cells_use,
  group.by = "predicted.id",
  reduction = "umap",
  label = FALSE,
  repel = TRUE,
  pt.size = 0.01,
  cols = My_Panel2
)
PPRDICT

# 07. Visualize label distributions and chromatin coverage ---------------------

# 07.1. Calculate ATAC-cluster proportions within each selected predicted type.
celltype_use <- names(My_Panel2)
heat_df <- ATAC@meta.data %>%
  dplyr::filter(predicted.id %in% celltype_use) %>%
  dplyr::count(predicted.id, seurat_clusters) %>%
  dplyr::group_by(predicted.id) %>%
  dplyr::mutate(prop = n / sum(n))
heat_df$predicted.id <- factor(
  heat_df$predicted.id,
  levels = c(
    "GZMK+ memory T", "memory T(2)", "memory T(1)", "exhausted T(1)"
  )
)

ggplot(heat_df, aes(seurat_clusters, predicted.id, fill = prop)) +
  geom_tile() +
  scale_fill_gradient(low = "grey95", high = "#2166AC") +
  theme_classic()

# 07.2. Plot coverage at selected genes and an mm10 genomic interval.
Idents(ATAC) <- "predicted.id"
# Set the plotting order on a subset to preserve all labels in the full object.
ATAC_coverage <- subset(ATAC, cells = cells_use)
DefaultAssay(ATAC_coverage) <- "peaks"
ATAC_coverage$predicted.id <- factor(
  ATAC_coverage$predicted.id,
  levels = c(
    "exhausted T(1)", "memory T(1)", "memory T(2)", "GZMK+ memory T"
  )
)
CoveragePlot(
  ATAC_coverage,
  cells = cells_use,
  "Gzmb",
  extend.downstream = 1000,
  group.by = "predicted.id"
) &
  scale_fill_manual(values = My_Panel2)
CoveragePlot(
  ATAC_coverage,
  cells = cells_use,
  "Gzma",
  extend.downstream = 2000,
  group.by = "predicted.id"
) &
  scale_fill_manual(values = My_Panel2)
CoveragePlot(
  ATAC_coverage,
  cells = cells_use,
  "chr13-113171630-113184500",
  group.by = "predicted.id"
) &
  scale_fill_manual(values = My_Panel2)

# 08. Add JASPAR motifs and calculate chromVAR deviations ----------------------

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(
    collection = "CORE",
    tax_group = "vertebrates",
    all_versions = FALSE
  )
)
DefaultAssay(ATAC) <- "peaks"
ATAC <- AddMotifs(
  object = ATAC,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

# Retain peaks on standard chromosomes and exclude mitochondrial peaks.
gr <- granges(ATAC[["peaks"]])
gr_std <- keepStandardChromosomes(gr, pruning.mode = "coarse")
gr_std <- gr_std[as.character(seqnames(gr_std)) != "chrM"]
keep_chr <- unique(as.character(seqnames(gr_std)))
keep_peaks <- rownames(ATAC[["peaks"]])[
  as.character(seqnames(gr)) %in% keep_chr
]
ATAC <- subset(ATAC, features = keep_peaks)

ATAC <- RunChromVAR(
  object = ATAC,
  genome = BSgenome.Mmusculus.UCSC.mm10
)
DefaultAssay(ATAC) <- "chromvar"

# 09. Identify differential motif activity in selected cell types --------------

# The original motif comparison includes GZMK+ inflammatory T, whereas the
# coverage plots in Section 07 use GZMK+ memory T. These selections are distinct.
Idents(ATAC) <- "predicted.id"
cells_use <- WhichCells(
  ATAC,
  idents = c(
    "exhausted T(1)", "memory T(1)", "memory T(2)", "GZMK+ inflammatory T"
  )
)
ATAC_sub <- subset(ATAC, cells = cells_use)
ATAC_sub$predicted.id <- factor(
  ATAC_sub$predicted.id,
  levels = c(
    "exhausted T(1)", "memory T(1)", "memory T(2)", "GZMK+ inflammatory T"
  )
)

motif_markers_all <- FindAllMarkers(
  ATAC_sub,
  only.pos = TRUE,
  logfc.threshold = 0,
  min.pct = 0.05,
  test.use = "wilcox"
)
top_motifs_per_cluster <- motif_markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()

heatmap_features <- unique(top_motifs_per_cluster$gene)
length(heatmap_features)

# Retain the original AverageExpression summary for reproducibility; review
# its suitability for chromVAR deviations separately before reusing this step.
avg_motif <- AverageExpression(
  ATAC_sub,
  assays = "chromvar",
  features = heatmap_features,
  group.by = "predicted.id",
  slot = "data"
)$chromvar
avg_motif_z <- t(scale(t(avg_motif)))
avg_motif_z[is.na(avg_motif_z)] <- 0

# 10. Annotate motif IDs and visualize standardized motif summaries ------------

pfm_list <- getMatrixSet(
  JASPAR2020,
  opts = list(collection = "CORE", matrixtype = "PFM", all_versions = TRUE)
)
anno_df <- tibble(
  motif_id = names(pfm_list),
  tf_name = sapply(pfm_list, name),
  tf_class = sapply(pfm_list, function(x) x@tags$class[1] %||% NA),
  tf_family = sapply(pfm_list, function(x) x@tags$family[1] %||% NA)
)

motif_ids <- rownames(avg_motif_z)
motif2tf <- anno_df %>%
  distinct(motif_id, .keep_all = TRUE) %>%
  dplyr::select(motif_id, tf_name)
tf_names <- motif2tf$tf_name[match(motif_ids, motif2tf$motif_id)]
rownames(avg_motif_z) <- tf_names

pheatmap(
  t(avg_motif_z),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  fontsize = 10,
  angle_col = 45,
  border_color = "black"
)

# 11. Identify differentially accessible peaks in ATAC cluster 1 ---------------

DefaultAssay(ATAC) <- "peaks"
Idents(ATAC) <- "seurat_clusters"

# Compare ATAC cluster 1 with all other cells using logistic regression.
cluster1.da <- FindMarkers(
  object = ATAC,
  ident.1 = "1",
  ident.2 = NULL,
  only.pos = FALSE,
  test.use = "LR",
  min.pct = 0.05,
  logfc.threshold = 0.5
)
head(cluster1.da)

# 11.1. Select peaks with increased accessibility using study-specific thresholds.
cluster1.open.peaks <- cluster1.da %>%
  rownames_to_column("peak") %>%
  mutate(pct_difference = pct.1 - pct.2) %>%
  dplyr::filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.5,
    pct.1 >= 0.20,
    pct_difference >= 0.10
  ) %>%
  dplyr::arrange(desc(avg_log2FC), desc(pct_difference))

# 11.2. Select peaks with decreased accessibility using study-specific thresholds.
cluster1.closed.peaks <- cluster1.da %>%
  rownames_to_column("peak") %>%
  mutate(pct_difference = pct.1 - pct.2) %>%
  dplyr::filter(
    p_val_adj < 0.05,
    avg_log2FC < -0.5,
    pct.2 >= 0.20,
    pct_difference <= -0.10
  ) %>%
  dplyr::arrange(avg_log2FC, pct_difference)

cluster1.open.peaks$status <- "open"
cluster1.closed.peaks$status <- "closed"
cluster1.diff.peaks <- dplyr::bind_rows(
  cluster1.open.peaks,
  cluster1.closed.peaks
)

# 12. Annotate differential peaks and summarize associated genes ---------------

cluster1.peak.annotation <- ClosestFeature(
  object = ATAC,
  regions = cluster1.diff.peaks$peak
) %>%
  as.data.frame()

# Join annotations to both increased- and decreased-accessibility peaks.
cluster1.peak.annotation <- cluster1.diff.peaks %>%
  left_join(cluster1.peak.annotation, by = c("peak" = "query_region")) %>%
  arrange(desc(avg_log2FC), p_val_adj)

# Preserve the original distance < 4,000 bp criterion without substituting
# abs(distance); the effect depends on the annotation's distance convention.
cluster1.peak.annotation <- cluster1.peak.annotation[
  !is.na(cluster1.peak.annotation$distance) &
    cluster1.peak.annotation$distance < 4000,
]

# Summarize the retained differential peaks by associated gene.
cluster1.peak.gene.group <- cluster1.peak.annotation %>%
  dplyr::filter(!is.na(gene_name)) %>%
  mutate(abs_distance = abs(distance)) %>%
  group_by(gene_name) %>%
  summarise(
    n_specific_peaks = n_distinct(peak),
    best_peak = peak[which.max(avg_log2FC)],
    max_peak_log2FC = max(avg_log2FC, na.rm = TRUE),
    max_pct_difference = max(pct_difference, na.rm = TRUE),
    min_abs_distance = min(abs_distance, na.rm = TRUE),
    min_peak_padj = min(p_val_adj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    desc(n_specific_peaks),
    desc(max_peak_log2FC),
    desc(max_pct_difference)
  )

# Select the 25 highest and 25 lowest gene-level maximum peak log2 fold changes.
# The legacy group strings are retained; they are not used in the final plot.
top25_genes <- cluster1.peak.gene.group %>%
  slice_max(order_by = max_peak_log2FC, n = 25, with_ties = FALSE) %>%
  mutate(group = "Top20_open")
bottom25_genes <- cluster1.peak.gene.group %>%
  slice_min(order_by = max_peak_log2FC, n = 25, with_ties = FALSE) %>%
  mutate(group = "Bottom20_closed")

heat_df <- bind_rows(top25_genes, bottom25_genes) %>%
  arrange(desc(max_peak_log2FC)) %>%
  dplyr::select(gene_name, max_peak_log2FC) %>%
  column_to_rownames("gene_name")

# 13. Compare peak-associated genes with differential RNA expression -----------

# Compare RNA cluster 6 with the other cells in the selected RNA reference.
Idents(RNA) <- "seurat_clusters"
cluster6.da <- FindMarkers(
  object = RNA,
  ident.1 = "6",
  ident.2 = NULL,
  only.pos = FALSE,
  test.use = "LR",
  min.pct = 0.05,
  logfc.threshold = 0
)

atac_df <- heat_df %>%
  as.data.frame() %>%
  rownames_to_column("gene_name")
rna_df <- cluster6.da %>%
  as.data.frame() %>%
  rownames_to_column("gene_name") %>%
  dplyr::select(
    gene_name,
    RNA_log2FC = avg_log2FC,
    RNA_padj = p_val_adj,
    RNA_pct1 = pct.1,
    RNA_pct2 = pct.2
  )
plot_df <- atac_df %>%
  left_join(rna_df, by = "gene_name")
plot_df <- plot_df %>%
  dplyr::mutate(
    RNA_direction = dplyr::case_when(
      RNA_log2FC > 0 ~ "Up",
      RNA_log2FC < 0 ~ "Down",
      TRUE ~ NA_character_
    ),
    RNA_significance = dplyr::case_when(
      is.na(RNA_log2FC) ~ NA_character_,
      !is.na(RNA_padj) & RNA_padj < 0.05 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )
plot_df$gene_name <- factor(
  plot_df$gene_name,
  levels = rev(rownames(heat_df))
)

# 14. Plot accessibility and expression changes side by side -------------------

# 14.1. Display the maximum peak log2 fold change associated with each gene.
p_atac <- ggplot(
  plot_df,
  aes(x = max_peak_log2FC, y = gene_name, fill = max_peak_log2FC > 0)
) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(
    values = c("TRUE" = "#D6A3A3", "FALSE" = "#91B3C4"),
    guide = "none"
  ) +
  labs(x = "ATAC peak log2FC", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 10)
  )

# 14.2. Display RNA fold changes, significance, and detection fractions.
gene_order <- rownames(heat_df)
p_rna <- ggplot(plot_df, aes(x = RNA_log2FC, y = gene_name)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +

  # Segment color indicates the direction of the RNA expression change.
  geom_segment(
    data = plot_df %>% dplyr::filter(!is.na(RNA_log2FC)),
    aes(
      x = 0,
      xend = RNA_log2FC,
      y = gene_name,
      yend = gene_name,
      color = RNA_direction
    ),
    linewidth = 0.7
  ) +

  # Filled and open circles distinguish significant and nonsignificant results.
  geom_point(
    data = plot_df %>% dplyr::filter(!is.na(RNA_log2FC)),
    aes(color = RNA_direction, shape = RNA_significance, size = RNA_pct1),
    stroke = 1
  ) +

  # Use a gray open circle where no RNA differential-expression result is present.
  geom_point(
    data = plot_df %>% dplyr::filter(is.na(RNA_log2FC)),
    aes(x = 0, y = gene_name),
    shape = 1,
    size = 2,
    stroke = 0.8,
    color = "grey60",
    inherit.aes = FALSE
  ) +
  scale_y_discrete(limits = rev(gene_order), drop = FALSE) +
  scale_color_manual(
    values = c("Up" = "#D6A3A3", "Down" = "#91B3C4"),
    name = "RNA direction"
  ) +
  scale_shape_manual(
    values = c("Significant" = 16, "Not significant" = 1),
    name = "RNA significance",
    labels = c(
      "Significant" = "Adjusted p < 0.05",
      "Not significant" = "Adjusted p ≥ 0.05"
    )
  ) +
  scale_size_continuous(range = c(1.8, 4), name = "RNA pct.1") +
  labs(x = "RNA expression log2FC", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 10),
    legend.position = "right"
  )
p_atac + p_rna
