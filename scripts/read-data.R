library(readr); library(dplyr); library(here)

# extract metadata
meta <- readr::read_tsv(here("data", "meta.tsv"))
glimpse(meta)
count(meta, biosample_id)
sce  <- zellkonverter::readH5AD(here("data", "anndata.h5ad"))

# compare nUMI and UMICount
colnames(meta)
meta$diff <- meta$nUMI - meta$`UMI Count`
summary(meta$diff)
mean(meta$diff == 0)                  # fraction of cells where they match
all(meta$`UMI Count` == round(meta$`UMI Count`))
