meta$diff <- meta$nUMI - meta$UMICount
summary(meta$diff)
mean(meta$diff == 0)                  # fraction of cells where they match
all(meta$UMICount == round(meta$UMICount))
min(meta$nUMI)
