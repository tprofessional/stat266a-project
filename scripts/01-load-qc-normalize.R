# Steps 1-3: load SHARE-seq RNA counts, QC, normalize.
# Output: outputs/sce_qc_norm.rds (SingleCellExperiment) + QC figures in figures/01-*.

library(here)
library(dplyr)
library(readr)
library(ggplot2)
library(Matrix)
library(SingleCellExperiment)
library(scuttle)
library(scran)
library(BiocParallel)

set.seed(266)
bpp <- MulticoreParam(workers = 8)
dir.create(here("outputs"), showWarnings = FALSE)
dir.create(here("figures"), showWarnings = FALSE)

# ---- 1. Load ---------------------------------------------------------------

# anndata.h5ad only holds 3,183 log-normalized HVGs, so read raw counts for all
# genes from the dense TSV (genes x cells), streaming it into a sparse matrix.
read_dense_tsv_sparse <- function(path, chunk_size = 2000) {
  con <- gzfile(path, "r")
  on.exit(close(con))
  cells <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]][-1]
  genes <- character(0)
  chunks <- list()
  repeat {
    lines <- readLines(con, n = chunk_size)
    if (length(lines) == 0) break
    dt <- data.table::fread(text = lines, sep = "\t", header = FALSE,
                            colClasses = c("character", rep("numeric", length(cells))))
    genes <- c(genes, dt[[1]])
    chunks[[length(chunks) + 1]] <- as(as.matrix(dt[, -1]), "CsparseMatrix")
    message(sprintf("  read %d genes", length(genes)))
  }
  m <- do.call(rbind, chunks)
  dimnames(m) <- list(make.unique(genes), cells)
  m
}

counts <- read_dense_tsv_sparse(here("data", "RNA", "exprMatrix.tsv.gz"))

# Everything downstream (scran size factors, pseudo-bulk edgeR) assumes raw counts.
stopifnot("X is not raw integer counts" = all(counts@x == round(counts@x)))

meta <- read_tsv(here("data", "RNA", "meta.tsv"), show_col_types = FALSE) |>
  dplyr::rename(celltype_fine = CellType)  # S4Vectors masks rename()
stopifnot(all(colnames(counts) %in% meta$cellId))
meta <- meta[match(colnames(counts), meta$cellId), ]

celltype_map <- read_csv(here("metadata", "celltype_map.csv"), show_col_types = FALSE)
stopifnot("unmapped cell types in metadata/celltype_map.csv" =
            all(meta$celltype_fine %in% celltype_map$celltype_fine))

# biosample_id looks like "Share-1m-4" / "Share-d23-1"
age_days <- c(d23 = 23, `1m` = 30, `2m` = 60, `3m` = 90)
meta <- meta |>
  left_join(select(celltype_map, celltype_fine, celltype, lineage), by = "celltype_fine") |>
  mutate(
    timepoint = sub("^Share-([^-]+)-.*$", "\\1", biosample_id),
    timepoint = factor(timepoint, levels = names(age_days)),
    age_days = unname(age_days[as.character(timepoint)]),
    organoid = biosample_id
  )

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = S4Vectors::DataFrame(meta, row.names = meta$cellId)
)
message(sprintf("Loaded %d genes x %d cells", nrow(sce), ncol(sce)))

# ---- 2. QC -----------------------------------------------------------------

is_mito <- grepl("^MT-", rownames(sce))
message(sprintf("%d mitochondrial genes", sum(is_mito)))
sce <- addPerCellQCMetrics(sce, subsets = list(Mito = is_mito), BPPARAM = bpp)

# The authors already filtered on library size (floor at ~1000 UMI), so only
# filter on mito. SHARE-seq profiles nuclei and mito % is mostly < 1%, so
# per-organoid MAD thresholds land at ~0.2% at 3m and discard healthy cells.
# A fixed cutoff only removes the real high-mito tail (mostly at 1m).
max_mito_pct <- 5
sce$discard <- sce$subsets_Mito_percent > max_mito_pct
print(table(organoid = sce$organoid, discard = sce$discard))
print(table(celltype = sce$celltype, discard = sce$discard))

qc_long <- as.data.frame(colData(sce)) |>
  select(organoid, timepoint, discard, sum, detected, any_of("subsets_Mito_percent")) |>
  tidyr::pivot_longer(c(sum, detected, any_of("subsets_Mito_percent")),
                      names_to = "metric", values_to = "value")
p_qc <- ggplot(qc_long, aes(organoid, value, colour = discard)) +
  geom_jitter(size = 0.2, width = 0.3, alpha = 0.4) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_y_log10() +
  geom_hline(data = data.frame(metric = "subsets_Mito_percent", y = max_mito_pct),
             aes(yintercept = y), linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here("figures", "01-qc-metrics.png"), p_qc, width = 9, height = 9, dpi = 150)

sce <- sce[, !sce$discard]

# Drop genes seen in too few cells to model.
keep_genes <- rowSums(counts(sce) > 0) >= 10
sce <- sce[keep_genes, ]
message(sprintf("After QC: %d genes x %d cells", nrow(sce), ncol(sce)))

# ---- 3. Normalize ----------------------------------------------------------

# Pooled size factors (Lun et al. 2016); cluster within organoid so pools stay
# within one sample.
clusters <- quickCluster(sce, block = sce$organoid, BPPARAM = bpp)
sce <- computeSumFactors(sce, clusters = clusters, BPPARAM = bpp)
stopifnot(all(sizeFactors(sce) > 0))
sce <- logNormCounts(sce)

p_sf <- ggplot(data.frame(lib = colSums(counts(sce)) / mean(colSums(counts(sce))),
                          sf = sizeFactors(sce), timepoint = sce$timepoint),
               aes(lib, sf, colour = timepoint)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_abline(linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "library size factor", y = "deconvolution size factor") +
  theme_bw()
ggsave(here("figures", "01-size-factors.png"), p_sf, width = 6, height = 5, dpi = 150)

saveRDS(sce, here("outputs", "sce_qc_norm.rds"))
message("Saved outputs/sce_qc_norm.rds")

# ---- Summaries for later steps ---------------------------------------------

cell_counts <- as.data.frame(colData(sce)) |>
  dplyr::count(timepoint, organoid, celltype) |>  # Bioconductor masks count()
  group_by(organoid) |>
  mutate(prop = n / sum(n)) |>
  ungroup()
write_csv(cell_counts, here("outputs", "celltype_counts_per_organoid.csv"))

p_comp <- ggplot(cell_counts, aes(organoid, prop, fill = celltype)) +
  geom_col() +
  facet_grid(~timepoint, scales = "free_x", space = "free_x") +
  labs(y = "fraction of cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here("figures", "01-composition.png"), p_comp, width = 10, height = 5, dpi = 150)
