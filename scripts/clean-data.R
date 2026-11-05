library(dplyr); library(stringr); library(here)
library(SummarizedExperiment); library(scuttle)
count <- dplyr::count

# 1. Load data and check size/normalization
meta <- readr::read_tsv(here("data", "RNA", "meta.tsv"), show_col_types = FALSE)
sce  <- zellkonverter::readH5AD(here("data", "RNA", "anndata.h5ad"))
dim(sce)                    # genes x cells
range(assay(sce))           # decimals = already log-normalized

# 2. Make duplicate gene names unique (GENE, GENE.1) so nothing is dropped
rownames(sce) <- make.unique(rownames(sce))

# 3. Line up metadata with cells (same order as the expression matrix)
meta <- meta[match(colnames(sce), meta$cellId), ]

# 4. Split biosample_id ("Share-1m-4") into timepoint, age, organoid
meta <- meta %>%
  mutate(
    timepoint = str_extract(biosample_id, "(d\\d+|\\d+m)"),
    organoid  = str_extract(biosample_id, "\\d+$"),
    age_days  = case_when(
      str_starts(timepoint, "d") ~ as.numeric(str_remove(timepoint, "d")),
      str_ends(timepoint, "m")   ~ as.numeric(str_remove(timepoint, "m")) * 30
    ),
    timepoint = factor(timepoint, levels = c("d23", "1m", "2m", "3m"))
  )
count(meta, biosample_id, timepoint, age_days, organoid)   # no NAs now

# 5. Standardize cell type labels
meta <- meta %>%
  mutate(cell_class = case_when(
    CellType %in% c("aRG", "oRG")                                  ~ "Radial glia",
    CellType %in% c("IP", "IPNewbornPN")                           ~ "Intermediate progenitors",
    CellType %in% c("Newborn CPN", "CPN", "CFuPN", "Immature DL PN",
                    "Unspecified PN", "Excitatory neu")            ~ "Excitatory neurons",
    CellType %in% c("Subcortical neurons", "Subcortical progenitors") ~ "Subcortical",
    CellType == "Cajal Retzius"                                    ~ "Cajal-Retzius",
    CellType %in% c("Cortical hem", "Mesenchyme", "Neural placode",
                    "PreplateSubplate")                            ~ "Other",
    CellType == "Unknown"                                          ~ NA_character_
  ))

count(meta, cell_class)                  # nothing unexpected should be NA besides Unknown
meta_clean <- meta %>% filter(!is.na(cell_class))
nrow(meta) - nrow(meta_clean)            # cells dropped as Unknown (report this)

# 6. Quality metrics: how many genes are detected per cell / cells per gene
cell_qc <- perCellQCMetrics(sce, assay.type = 1)
gene_qc <- perFeatureQCMetrics(sce, assay.type = 1)
meta$n_genes_detected <- cell_qc$detected
summary(gene_qc$detected)   # % of cells expressing each gene

# 7. Keep genes expressed (> 0) in at least 10 cells
keep <- Matrix::rowSums(assay(sce) > 0) >= 10
c(before = nrow(sce), after = sum(keep))
sce <- sce[keep, ]

summary(sce)
dim(sce)

# 8. Drop Unknown cells for the composition analysis
meta_clean <- meta %>% filter(!is.na(cell_class))
nrow(meta) - nrow(meta_clean)            # cells dropped as Unknown (report this)


comp <- meta_clean %>%
  count(biosample_id, timepoint, age_days, cell_class) %>%
  group_by(biosample_id) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

comp %>% group_by(biosample_id) %>% summarise(total = sum(prop))   # every row should be 1
head(comp)

library(ggplot2)
p_comp <- ggplot(comp, aes(x = biosample_id, y = prop, fill = cell_class)) +
  geom_col() +
  facet_grid(~ timepoint, scales = "free_x", space = "free_x") +
  labs(title = "Cell type composition shifts across organoid age",
       x = "Organoid (biosample)", y = "Proportion of cells",
       fill = "Cell class",
       caption = "Each bar is one organoid; panels group organoids by timepoint.") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here("figures", "clean-data-composition.png"), p_comp, width = 9, height = 5, dpi = 200)

  # Prepping for heatmap of top 40 genes

library(tidyr); library(limma) #tidyr for data cleaning and limma for gene expression analyis

# getting cells that we will use from clean set
cells_use <- meta_clean$cellId
# putting into a single cell experiment matrix
sce_use <- sce[, cells_use]
# metadata matching to cells being used
md_use  <- meta_clean[match(cells_use, meta_clean$cellId), ]
# calculating mean of gene expression per organoid
pb <- summarizeAssayByGroup(sce_use, ids = md_use$biosample_id,
                            statistics = "mean", assay.type = 1)
sample_info <- md_use %>%
  distinct(biosample_id, timepoint, age_days, organoid) %>%
  arrange(timepoint, organoid)

pb_mat <- assay(pb, "mean")[, sample_info$biosample_id]   # genes x organoids, ordered by age

#Rank genes by change across timepoints (F-test across all timepoints)
expressed <- rowMeans(pb_mat) > 0.1          # dropping barely expressed genes
design <- model.matrix(~ timepoint, data = sample_info)
fit <- eBayes(lmFit(pb_mat[expressed, ], design))
top <- topTable(fit, coef = 2:ncol(design), number = Inf, sort.by = "F")
head(top)

# selecting the top 8 genes expressed
line_genes <- rownames(top)[1:8]

# reshaping the matrix into the long table format that ggplot is going to use
line_df <- as.data.frame(pb_mat[line_genes, , drop = FALSE]) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "biosample_id", values_to = "expr") %>%
  left_join(sample_info, by = "biosample_id") %>%
  mutate(gene = factor(gene, levels = line_genes))

#plotting timepoint on the x axis and expression on the y axis
p_genes <- ggplot(line_df, aes(timepoint, expr)) +
  geom_point(alpha = 0.6, position = position_jitter(width = 0.1, height = 0)) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 0.8) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  labs(title = "Change in expression of top-ranked genes over time",
       x = "Timepoint", y = "Mean expression per organoid (log-based)",
       caption = "Points are organoids, line connects timepoint means.") +
  theme_minimal()
ggsave(here("figures", "top8-genes-over-time.png"), p_genes, width = 10, height = 5, dpi = 200)

