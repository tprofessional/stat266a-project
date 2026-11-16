# run clean_data.R before this
library(tidyr)

# only keep cells with known types (meta_clean gets rid of Unknown type)
sce_clean <- sce[, meta_clean$cellId]

# for each gene compute the average in each organoid
cells_by_organoid <- split(seq_len(ncol(sce_clean)), meta_clean$biosample_id)
gene_means <- sapply(cells_by_organoid, function(i) Matrix::rowMeans(assay(sce_clean)[, i]))

# 
y_tbl <- comp %>%
  # add a row of 0s if an organoid has no cells in the class (prevents from losing a datapoint in early stage when there may not be any excitatory neurons)
  complete(biosample_id, cell_class, fill = list(n = 0, prop = 0)) %>%
  # pick the excitatory neurons class to model
  filter(cell_class == "Excitatory neurons") %>%
  # order the organoids to match the order of gene_means to ensure proportions are assigned to the correct organoid
  arrange(match(biosample_id, colnames(gene_means)))
# don't continue if the orders don't match
stopifnot(all(y_tbl$biosample_id == colnames(gene_means)))
# vector of 12 proportions
y <- y_tbl$prop

# head(y)

# perform regression per gene
fit_gene <- function(x) {
  x <- as.numeric(scale(x))
  s <- summary(lm(y ~ x))
  c(slope = s$coefficients["x", "Estimate"],
    p     = s$coefficients["x", "Pr(>|t|)"],
    r2    = s$r.squared)
}

results <- t(apply(gene_means, 1, fit_gene)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(desc(r2))

# top 10 genes
head(results, 10)

# top gene
top_gene <- results$gene[1]
model_df <- data.frame(prop = y, expr = gene_means[top_gene, ],
                       timepoint = y_tbl$timepoint)
m <- lm(prop ~ expr, data = model_df)

summary(m)
confint(m)

# residual diagnostics
plot(m, which = 1)
shapiro.test(resid(m))
