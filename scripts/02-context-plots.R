# Step 4 + context figures: UMAP, composition bars, cell type proportions over age.
# Input:  outputs/sce_qc_norm.rds (from 01-load-qc-normalize.R)
# Output: outputs/sce_dimred.rds, outputs/composition_trend_test.csv, figures/02-*

library(here)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(SingleCellExperiment)
library(scran)
library(scater)
library(limma)
library(BiocParallel)

set.seed(266)
bpp <- MulticoreParam(workers = 8)

sce <- readRDS(here("outputs", "sce_qc_norm.rds"))

# Fixed order and colors so every figure in the project matches.
celltype_levels <- c("aRG", "oRG", "IP", "DL_PN", "UL_PN", "Unspecified_PN",
                     "Cajal_Retzius", "Cortical_hem", "Non_cortical", "Unknown")
celltype_colors <- c(
  aRG = "#1b9e77", oRG = "#66c2a5", IP = "#a6d854",
  DL_PN = "#3a5fcd", UL_PN = "#8da0cb", Unspecified_PN = "#e78ac3",
  Cajal_Retzius = "#e6ab02", Cortical_hem = "#a6761d",
  Non_cortical = "#fc8d62", Unknown = "grey70"
)
timepoint_colors <- c(d23 = "#fdcc8a", `1m` = "#fc8d59", `2m` = "#e34a33", `3m` = "#7f0000")

sce$celltype <- factor(sce$celltype, levels = celltype_levels)
# Size factors split d23+3m vs 1m+2m (figures/01-size-factors.png), which points
# to two sequencing batches. Label it so it can be checked on every plot.
sce$seq_batch <- factor(ifelse(sce$timepoint %in% c("d23", "3m"), "d23 + 3m", "1m + 2m"))

# ---- 4. Dimensionality reduction -------------------------------------------

# PCA + UMAP take ~10 min; reuse them unless step 01 has been rerun since.
dimred_path <- here("outputs", "sce_dimred.rds")
if (file.exists(dimred_path) &&
    file.mtime(dimred_path) > file.mtime(here("outputs", "sce_qc_norm.rds"))) {
  sce_dimred <- readRDS(dimred_path)
  reducedDims(sce) <- reducedDims(sce_dimred)[, colnames(sce)]
  rm(sce_dimred)
} else {
  gene_var <- modelGeneVar(sce, block = sce$organoid, BPPARAM = bpp)
  hvgs <- getTopHVGs(gene_var, n = 2000)
  hvgs <- hvgs[!grepl("^MT-", hvgs)]

  sce <- runPCA(sce, ncomponents = 30, subset_row = hvgs, BPPARAM = bpp)
  sce <- runUMAP(sce, dimred = "PCA", n_neighbors = 30, BPPARAM = bpp)
  saveRDS(sce, dimred_path)
}

umap <- reducedDim(sce, "UMAP")
umap_df <- data.frame(UMAP1 = umap[, 1], UMAP2 = umap[, 2], as.data.frame(colData(sce))) |>
  slice_sample(prop = 1)  # shuffle so no group is always drawn on top

umap_plot <- function(colour_by, colours, title) {
  ggplot(umap_df, aes(UMAP1, UMAP2, colour = .data[[colour_by]])) +
    geom_point(size = 0.15, alpha = 0.6) +
    scale_colour_manual(values = colours, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = title) +
    coord_equal() +
    theme_void() +
    theme(plot.title = element_text(face = "bold"))
}

p_umap <- umap_plot("celltype", celltype_colors, "Cell type") +
  umap_plot("timepoint", timepoint_colors, "Organoid age") +
  umap_plot("seq_batch", c(`d23 + 3m` = "#4d9221", `1m + 2m` = "#c51b7d"),
            "Sequencing batch (inferred)") +
  plot_layout(ncol = 3)
ggsave(here("figures", "02-umap.png"), p_umap, width = 16, height = 5.5, dpi = 200)

# ---- Composition per organoid ----------------------------------------------

comp <- as.data.frame(colData(sce)) |>
  dplyr::count(timepoint, age_days, organoid, celltype, .drop = FALSE) |>
  group_by(organoid) |>
  filter(sum(n) > 0) |>
  mutate(prop = n / sum(n)) |>
  ungroup()
# .drop = FALSE crosses every organoid with every timepoint; keep real pairs
comp <- comp |> semi_join(distinct(as.data.frame(colData(sce)), organoid, timepoint),
                          by = c("organoid", "timepoint"))

p_bars <- ggplot(comp, aes(organoid, prop, fill = celltype)) +
  geom_col(width = 0.85) +
  facet_grid(~timepoint, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = celltype_colors, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_x_discrete(labels = \(x) sub("^Share-[^-]+-", "org ", x)) +
  labs(x = NULL, y = "Cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95"))
ggsave(here("figures", "02-composition-bars.png"), p_bars, width = 9, height = 4.5, dpi = 200)

# ---- Trend test: does each cell type's proportion change with age? ----------

# propeller-style (Phipson et al. 2022): arcsine-sqrt transformed proportions,
# one limma linear model per cell type with organoid as the replicate.
prop_mat <- comp |>
  select(organoid, celltype, prop) |>
  tidyr::pivot_wider(names_from = organoid, values_from = prop) |>
  tibble::column_to_rownames("celltype") |>
  as.matrix()
organoid_info <- distinct(comp, organoid, timepoint, age_days)
organoid_info <- organoid_info[match(colnames(prop_mat), organoid_info$organoid), ]
y <- asin(sqrt(prop_mat))

# Any change across timepoints (F-test). Catches rise-then-fall patterns
# (oRG peaks at 2m, IP at 1m) that a straight line misses.
design_tp <- model.matrix(~ timepoint, data = organoid_info)
fit_tp <- eBayes(lmFit(y, design_tp), robust = TRUE)
any_change <- topTable(fit_tp, coef = 2:ncol(design_tp), number = Inf, sort.by = "none") |>
  tibble::rownames_to_column("celltype") |>
  select(celltype, F, FDR_any_change = adj.P.Val)

# Monotonic direction: linear slope in age.
design_lin <- model.matrix(~ age_days, data = organoid_info)
fit_lin <- eBayes(lmFit(y, design_lin), robust = TRUE)
linear <- topTable(fit_lin, coef = "age_days", number = Inf, sort.by = "none") |>
  tibble::rownames_to_column("celltype") |>
  select(celltype, slope_per_day = logFC, FDR_linear = adj.P.Val)

trend <- left_join(any_change, linear, by = "celltype")
write_csv(trend, here("outputs", "composition_trend_test.csv"))
print(trend)

stars <- function(p) cut(p, c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "ns"))
facet_labels <- trend |>
  mutate(label = sprintf("%s  (%s)", celltype, stars(FDR_any_change))) |>
  select(celltype, label) |>
  tibble::deframe()

p_trend <- ggplot(comp, aes(age_days, prop, colour = celltype)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  geom_point(size = 1.8, alpha = 0.8) +
  facet_wrap(~celltype, scales = "free_y", ncol = 5,
             labeller = as_labeller(facet_labels)) +
  scale_colour_manual(values = celltype_colors, guide = "none") +
  scale_x_continuous(breaks = c(23, 30, 60, 90), labels = c("d23", "1m", "2m", "3m"),
                     guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Organoid age", y = "Fraction of cells",
       caption = "Points = organoids; line = timepoint mean. Stars: FDR for any change across timepoints (F-test on arcsine-sqrt proportions, limma).") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey95"),
        panel.grid.minor = element_blank())
ggsave(here("figures", "02-proportions-over-age.png"), p_trend, width = 13, height = 5.5, dpi = 200)
