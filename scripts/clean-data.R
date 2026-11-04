library(dplyr); library(SummarizedExperiment)

# align meta to sce
all(colnames(sce) %in% meta$cellId)    # want TRUE
anyDuplicated(meta$cellId)             # want 0
meta <- meta[match(colnames(sce), meta$cellId), ]
all(colnames(sce) == meta$cellId)      # want TRUE

# row counts
dim(sce); nrow(meta)
count(meta, biosample_id)
count(meta, CellType, sort = TRUE)

# missing values
colSums(is.na(meta))

# oddities
assayNames(sce); range(assay(sce))     # confirms normalized, not raw counts
sum(rowSums(assay(sce)) == 0)          # genes with zero expression
summary(meta$nUMI); summary(meta$nGene)