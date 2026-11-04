library(dplyr); library(stringr); library(here)
library(SummarizedExperiment); library(scuttle)
count <- dplyr::count

# 1. Load data and check size/normalization
meta <- readr::read_tsv(here("data", "meta.tsv"))
sce  <- zellkonverter::readH5AD(here("data", "anndata.h5ad"))
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
ggplot(comp, aes(x = biosample_id, y = prop, fill = cell_class)) +
  geom_col() +
  facet_grid(~ timepoint, scales = "free_x", space = "free_x") +
  labs(title = "Cell type composition shifts across organoid age",
       x = "Organoid (biosample)", y = "Proportion of cells",
       fill = "Cell class",
       caption = "Each bar is one organoid; panels group organoids by timepoint.") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))