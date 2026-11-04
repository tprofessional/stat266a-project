library(dplyr); library(tidyr); library(stringr); library(here)
library(SummarizedExperiment)

# 1. Load + check dims / normalization
meta <- readr::read_tsv(here("data", "meta.tsv"))
sce  <- zellkonverter::readH5AD(here("data", "anndata.h5ad"))
dim(sce)                      # genes x cells
range(assay(sce))             # small non-integer values = log-normalized

# 2. Make gene names unique (GENE, GENE.1, ...) without dropping any rows
sum(duplicated(rownames(sce)))
rownames(sce) <- make.unique(rownames(sce))

# 3. Attach metadata by cellId
stopifnot(all(colnames(sce) %in% meta$cellId), !anyDuplicated(meta$cellId))
meta <- meta[match(colnames(sce), meta$cellId), ]
stopifnot(all(colnames(sce) == meta$cellId))
colData(sce) <- cbind(colData(sce), DataFrame(meta[, -1]))   # cellId is already the colnames

head(meta)
dim(sce)

# 4. Parse biosample_id ("Share-1m-4") -> timepoint, age_days, organoid
meta <- meta %>%
  mutate(
    timepoint = str_extract(biosample_id, "\\d+[dm]"),
    organoid  = str_extract(biosample_id, "\\d+$"),
    age_days  = as.numeric(str_extract(timepoint, "\\d+")) *
                if_else(str_detect(timepoint, "m$"), 30, 1)
  )
dplyr::count(meta, biosample_id, timepoint, age_days, organoid)   # eyeball this

# 5. Standardize cell type labels (see options below)
dplyr::count(meta, timepoint, CellType) %>% pivot_wider(names_from = timepoint, values_from = n)

# 6. Per-cell QC
meta <- meta %>%
  mutate(
    nGene_check = colSums(assay(sce) > 0),
    mito_expr   = colMeans(assay(sce)[str_detect(rownames(sce), "^MT-"), , drop = FALSE]),
    log10_nUMI  = log10(nUMI)
  ) %>%
  select(-`UMI Count`, -Cluster, -`Louvain Cluster`)
summary(meta[, c("nGene", "nGene_check", "nUMI", "mito_expr")])
