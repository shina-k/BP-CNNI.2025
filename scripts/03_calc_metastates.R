library(tidyverse)
library(data.table)

ave_data <- fread("your_datasets_occupancy.csv")
ave_d1 <- fread("your_datasets_occupancy_ms1.csv") 
ave_d2 <- fread("your_datasets_occupancy_ms2.csv") 
tmat_group1 <- fread("your_datasets_transition_ms1.csv")
tmat_group2 <- fread("your_datasets_transition_ms2.csv")

################################################################################
#Meta-state estimation
#PCA
pca <- prcomp(ave_data)

loadings <- as.data.frame(pca$rotation[, 1])

# Hierarchical clustering
dist_matrix <- dist(loadings)  
hc <- hclust(dist_matrix, method = "ward.D2") 
plot(hc, main = "Dendrogram using Ward's method", xlab = "", sub = "")
dend <- hc %>% as.dendrogram()

cluster <- cutree(hc,k=3)
n_cluster <- length(table(cluster))
cluster <- cluster[order.dendrogram(dend)]

################################################################################
#Calc Entropy
calculate_entropy <- function(tmat, avd) {
  n <- sqrt(colnames(tmat) %>% length())
  transition_matrix <- matrix(as.numeric(tmat), nrow = n, byrow = TRUE)
  transition_matrix <- transition_matrix / rowSums(transition_matrix)
  d <- avd / sum(avd)
  row_entropy <- apply(transition_matrix, 1, function(row) {
    row <- row[row > 0]  
    -sum(row * log(row))
  })
  total_entropy <- sum(d * row_entropy)
  return(total_entropy)
}

x <- data.frame(x = 1:nrow(tmat_group1) )
entropy_results1 <- apply(x, 1, function(i) {
  calculate_entropy(tmat_group1[i, -1], ave_d1[i, -1])
})

entropy_results2 <- apply(x, 1, function(i) {
  calculate_entropy(tmat_group2[i, -1], ave_d2[i, -1])
})

