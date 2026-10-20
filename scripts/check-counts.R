library(zellkonverter)
library(SingleCellExperiment)

sce <- readH5AD("data/anndata.h5ad")
assayNames(sce)                 # look for "counts" or a raw layer
x <- assay(sce, 1)[1:200, 1:200]
all(x == round(x))              # TRUE = integer counts, FALSE = normalized
