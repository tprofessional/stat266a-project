library(readr); library(dplyr); library(here)

meta <- readr::read_tsv(here("data", "RNA", "meta.tsv"))
glimpse(meta)
count(meta, biosample_id)
