# =============================================================================
# Mouse scRNA-seq and paired TCR V(D)J analysis
# =============================================================================
# Purpose:
#   Perform quality control, clustering, cell-type annotation, visualization,
#   and paired TCR clonotype assignment for four mouse samples.
#
# Inputs:
#   - 10x Genomics feature-barcode matrices containing "Gene Expression".
#   - Per-sample filtered_contig_annotations.csv files from V(D)J analysis.
#
# Project layout (all paths are relative to the project root):
#   data/gene_expression/Sample1/ ... Sample4/
#   data/vdj/Sample1/filtered_contig_annotations.csv ... Sample4/
#   results/
#
# Usage:
#   Run from the project root after configuring the paths and sample labels
#   below. Review cluster selections and cell-type labels for other datasets.
#   Use a Seurat version supporting JoinLayers() and the layer argument.
# =============================================================================

# 01. Load required packages ---------------------------------------------------

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(readr)
library(tibble)
library(scales)

# 02. Configure input paths, output paths, and sample labels --------------------

# Edit these relative paths to match the local project layout.
base_dir <- file.path("data", "gene_expression")
vdj_base <- file.path("data", "vdj")
output_dir <- "results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_map <- c(
  "Sample1" = "2dpi_Ai14pos",
  "Sample2" = "2dpi_Ai14neg",
  "Sample3" = "14dpi_Ai14pos",
  "Sample4" = "14dpi_Ai14neg"
)

# 03. Import gene-expression matrices ----------------------------------------

#' Read a sample matrix and create a Seurat object.
#'
#' @param sample_id Name of the sample subdirectory under base_dir.
#' @param sample_label Descriptive label stored in the sample metadata.
#' @return A Seurat object containing counts and sample identifiers.
read_one <- function(sample_id, sample_label) {
  mat_dir <- file.path(base_dir, sample_id)
  if (!dir.exists(mat_dir)) stop("Directory not found: ", mat_dir)

  x <- Read10X(data.dir = mat_dir)[["Gene Expression"]]
  obj <- CreateSeuratObject(
    counts = x,
    project = sample_label,
    min.cells = 3,
    min.features = 200
  )
  obj$sample <- sample_label
  obj$sample_id <- sample_id
  obj
}

obj_list <- mapply(
  read_one,
  names(sample_map),
  sample_map,
  SIMPLIFY = FALSE
)
sapply(obj_list, ncol)

# 04. Calculate quality-control metrics and filter cells -----------------------

# Identify mouse mitochondrial genes by the "mt-" prefix.
obj_list <- lapply(obj_list, function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj
})

lapply(names(obj_list), function(nm) {
  VlnPlot(
    obj_list[[nm]],
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    layer = "counts",
    ncol = 3
  ) +
    ggtitle(nm)
})

#' Apply the study-specific cell quality-control thresholds.
#'
#' @param obj A Seurat object with percent.mt in its metadata.
#' @return A Seurat object retaining cells with 300-4,000 detected genes,
#'   at least 500 RNA counts, and no more than 5% mitochondrial counts.
qc_filter_one <- function(obj) {
  subset(
    obj,
    subset =
      nFeature_RNA >= 300 &
      nFeature_RNA <= 4000 &
      nCount_RNA >= 500 &
      percent.mt <= 5
  )
}

obj_list_qc <- lapply(obj_list, qc_filter_one)
rbind(
  before = sapply(obj_list, ncol),
  after = sapply(obj_list_qc, ncol)
)

# 05. Merge samples and perform initial dimensionality reduction --------------

sample.merge <- merge(
  x = obj_list_qc[[1]],
  y = obj_list_qc[-1],
  add.cell.ids = names(sample_map),
  project = "Poly_BALF_CD8"
)
sample.merge <- NormalizeData(sample.merge)
sample.merge <- FindVariableFeatures(
  sample.merge,
  selection.method = "vst",
  nfeatures = 3000
)
sample.merge <- ScaleData(sample.merge)
sample.merge <- RunPCA(sample.merge, npcs = 30)
sample.merge <- RunUMAP(sample.merge, reduction = "pca", dims = 1:15)
DimPlot(sample.merge, split.by = "sample")

# 06. Identify initial clusters and select alpha-beta T cells -------------------

sample.merge <- FindNeighbors(sample.merge, dims = 1:15)
sample.merge <- FindClusters(sample.merge, resolution = 0.7)
sample.merge <- JoinLayers(sample.merge)

markers <- FindAllMarkers(
  sample.merge,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)
write.csv(
  top10,
  file = file.path(output_dir, "Poly_top10_cluster_markers.csv"),
  row.names = FALSE
)

# Cluster IDs are specific to the study dataset and require biological review.
Idents(sample.merge) <- "seurat_clusters"
ab_T_merge <- subset(sample.merge, idents = c(0, 1, 2, 3, 4, 5, 6, 7))

# 07. Exclude TCR-associated genes --------------------------------------------

# Remove genes matching the study's TCR pattern before reclustering.
all_genes <- rownames(ab_T_merge)
tcr_pat <- "(^Tr[abdg]v)|(^Tr[abdg]j)|(^Tr[abdg]c)|(^Trac$)|(^Trbc[12]?$)|(^Trdc$)|(^Tcr[abdg])"
selected_genes <- all_genes[
  !grepl(tcr_pat, all_genes, ignore.case = TRUE, perl = TRUE)
]
ab_T_merge_selected_genes <- subset(ab_T_merge, features = selected_genes)

# 08. Reanalyze the T-cell subset and export cluster markers -------------------

ab_T_merge_selected_genes <- NormalizeData(ab_T_merge_selected_genes)
ab_T_merge_selected_genes <- FindVariableFeatures(ab_T_merge_selected_genes)
ab_T_merge_selected_genes <- ScaleData(ab_T_merge_selected_genes)
ab_T_merge_selected_genes <- RunPCA(ab_T_merge_selected_genes)
ElbowPlot(ab_T_merge_selected_genes, ndims = 50)
ab_T_merge_selected_genes <- FindNeighbors(
  ab_T_merge_selected_genes,
  dims = 1:10
)
ab_T_merge_selected_genes <- FindClusters(
  ab_T_merge_selected_genes,
  resolution = 0.5
)
ab_T_merge_selected_genes <- RunUMAP(ab_T_merge_selected_genes, dims = 1:10)

# Save the clustered object before the optional cell-type annotation below.
saveRDS(
  ab_T_merge_selected_genes,
  file = file.path(output_dir, "ab_T_merge_selected_genes.rds")
)

markers <- FindAllMarkers(
  ab_T_merge_selected_genes,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)
write.csv(
  top10,
  file = file.path(output_dir, "Poly_CD8_markers.csv"),
  row.names = FALSE
)

# 09. Assign cell-type labels and visualize the T-cell subset -------------------

# Optional annotation: labels correspond to the current cluster-level order.
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
Idents(ab_T_merge_selected_genes) <- "seurat_clusters"
names(celltype) <- levels(ab_T_merge_selected_genes)
ab_T_merge_selected_genes <- RenameIdents(ab_T_merge_selected_genes, celltype)
ab_T_merge_selected_genes$celltype <- Idents(ab_T_merge_selected_genes)

# 09.1. Plot cell-type UMAP embeddings.
My_Panel <- c(
  "#a2d2e7", "#67abcd", "#3581b7", "#6fb3a8", "#b3e19b",
  "#50aa4b", "#ffc17f", "#cf9f88", "#ff8831"
)
names(My_Panel) <- levels(Idents(ab_T_merge_selected_genes))

emb <- Embeddings(ab_T_merge_selected_genes, "umap") %>%
  as.data.frame()
emb$cluster <- factor(Idents(ab_T_merge_selected_genes))
centers <- emb %>%
  group_by(cluster) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )

PDIM <- DimPlot(ab_T_merge_selected_genes, pt.size = 0.01, cols = My_Panel)
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
  ab_T_merge_selected_genes,
  pt.size = 0.01,
  cols = My_Panel,
  split.by = "sample"
) &
  theme(legend.position = "none")
PDIM_SEP +
  geom_point(
    data = centers,
    aes(umap_1, umap_2, colour = cluster),
    shape = 21,
    fill = alpha("white", 0.5),
    size = 6,
    stroke = 0.3
  ) +
  geom_text(
    data = centers,
    aes(umap_1, umap_2, label = cluster),
    color = "black",
    size = 4
  )

# 09.2. Calculate and plot cell-type proportions within each sample.
meta <- ab_T_merge_selected_genes@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(cluster = Idents(ab_T_merge_selected_genes))
df_prop <- meta %>%
  group_by(sample, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(freq = n / sum(n))

ggplot(df_prop, aes(x = sample, y = freq, fill = cluster)) +
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

# 09.3. Plot expression of selected CD8 T-cell markers.
markers_cd8 <- c(
  "Sell", "Ccr7", "Il7r", "Tcf7", "Lef1", "Eomes", "Klf2", "Bcl2",
  "Gzma", "Gzmb", "Gzmk", "Ly6a", "Epsti1", "Ccl5", "Ifng", "Prf1",
  "Fasl", "Pdcd1", "Ctla4", "Tox", "Tigit", "Lag3", "Havcr2",
  "Ccr2", "Cxcr3", "Cxcr6", "Itgae", "Itga4"
)
markers_Tbox <- c("Tbx2", "Eomes", "Tbx21", "Tbx20", "Tbr1")

p <- DotPlot(ab_T_merge_selected_genes, features = markers_cd8) +
  guides(color = "none")
p +
  geom_point(
    aes(size = pct.exp, fill = avg.exp.scaled),
    shape = 21,
    color = "black",
    stroke = 0.35
  ) +
  scale_fill_gradientn(
    colours = c("#3C4EC2", "#8DB0FE", "#F7F7F7", "#F4987A", "#B40426"),
    guide = guide_colorbar()
  ) +
  scale_size(
    range = c(0, 8),
    guide = guide_legend()
  ) +
  guides(color = "none") +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title = element_blank()
  )

# 10. Import and summarize TCR V(D)J annotations -------------------------------

sample_ids <- c("Sample1", "Sample2", "Sample3", "Sample4")

# Reload the clustered snapshot saved in Section 08, as in the original workflow.
# Optional cell-type labels assigned in Section 09 are absent from this snapshot.
ab_T_merge_selected_genes <- readRDS(
  file = file.path(output_dir, "ab_T_merge_selected_genes.rds")
)

#' Read high-confidence V(D)J contigs for one sample.
#'
#' @param sid Name of the sample subdirectory under vdj_base.
#' @return A contig-level data frame with the sample_id column added.
#' @details Contigs with missing full_length values are retained.
read_vdj_one <- function(sid) {
  f <- file.path(vdj_base, sid, "filtered_contig_annotations.csv")
  if (!file.exists(f)) stop("VDJ file not found: ", f)

  df <- readr::read_csv(f, show_col_types = FALSE)
  df <- df %>%
    filter(high_confidence == TRUE) %>%
    filter(is.na(full_length) | full_length == TRUE)

  df$sample_id <- sid
  df
}

vdj_raw_list <- lapply(sample_ids, read_vdj_one)
names(vdj_raw_list) <- sample_ids

#' Summarize V(D)J contigs into one record per sample and cell barcode.
#'
#' @param df A contig-level data frame containing Cell Ranger V(D)J annotations
#'   and a sample_id column.
#' @return A data frame with chain-pairing indicators, clonotype keys, CDR3
#'   amino acid and nucleotide sequences, and V(D)J gene assignments.
#' @details TRA keys use CDR3 nucleotide sequences and V/J genes. TRB keys use
#'   CDR3 nucleotide sequences and V/D/J genes. Multiple chain records are
#'   deduplicated, sorted, and concatenated before combining the TRA/TRB keys.
#'   Missing full_length and productive values are retained during filtering.
vdj_summarize_one <- function(df) {
  df %>%
    filter(is_cell == TRUE) %>%
    filter(high_confidence == TRUE) %>%
    filter(is.na(full_length) | full_length) %>%
    filter(is.na(productive) | productive) %>%
    group_by(sample_id, barcode) %>%
    summarise(
      vdj_n_contigs = n(),
      vdj_has_TRA = any(chain == "TRA"),
      vdj_has_TRB = any(chain == "TRB"),
      TRA_key = paste(
        sort(unique(na.omit(
          paste0(
            "cdr3nt=", cdr3_nt[chain == "TRA"],
            "|V=", v_gene[chain == "TRA"],
            "|J=", j_gene[chain == "TRA"]
          )
        ))),
        collapse = ";"
      ),
      TRB_key = paste(
        sort(unique(na.omit(
          paste0(
            "cdr3nt=", cdr3_nt[chain == "TRB"],
            "|V=", v_gene[chain == "TRB"],
            "|D=", d_gene[chain == "TRB"],
            "|J=", j_gene[chain == "TRB"]
          )
        ))),
        collapse = ";"
      ),
      clonotype_key = paste(
        paste0("TRA:", ifelse(TRA_key == "", "NA", TRA_key)),
        paste0("TRB:", ifelse(TRB_key == "", "NA", TRB_key)),
        sep = "|"
      ),
      raw_clonotype_id = dplyr::first(na.omit(raw_clonotype_id)),
      vdj_cdr3aa = paste(
        sort(unique(na.omit(paste0(chain, ":", cdr3)))),
        collapse = ";"
      ),
      vdj_cdr3nt = paste(
        sort(unique(na.omit(paste0(chain, ":", cdr3_nt)))),
        collapse = ";"
      ),
      vdj_vj = paste(
        sort(unique(na.omit(
          ifelse(
            chain == "TRA",
            paste0(chain, ":V=", v_gene, "|J=", j_gene),
            paste0(chain, ":V=", v_gene, "|D=", d_gene, "|J=", j_gene)
          )
        ))),
        collapse = ";"
      ),
      .groups = "drop"
    )
}

vdj_sum <- bind_rows(lapply(vdj_raw_list, vdj_summarize_one))

# 11. Assign clonotypes across samples and add Seurat metadata -----------------

# Require both TRA and TRB chains and assign a shared ID to each clonotype key.
# Clonotype numbering applies to the combined dataset analyzed in this run.
vdj_sum <- vdj_sum %>%
  filter(vdj_has_TRA & vdj_has_TRB) %>%
  arrange(clonotype_key) %>%
  group_by(clonotype_key) %>%
  mutate(clonotype_id = paste0("clonotype", cur_group_id())) %>%
  ungroup()

# Match the sample-prefixed cell names created during Seurat object merging.
vdj_sum2 <- vdj_sum %>%
  mutate(cell = paste0(sample_id, "_", barcode)) %>%
  select(
    cell, clonotype_id, vdj_n_contigs, vdj_has_TRA, vdj_has_TRB,
    vdj_cdr3aa, vdj_cdr3nt, vdj_vj, TRA_key, TRB_key, clonotype_key
  )

ab_T_merge_selected_genes <- AddMetaData(
  ab_T_merge_selected_genes,
  metadata = vdj_sum2 %>% tibble::column_to_rownames("cell")
)

# Report the numbers of cells with and without an assigned paired clonotype.
table(is.na(ab_T_merge_selected_genes$clonotype_id))
