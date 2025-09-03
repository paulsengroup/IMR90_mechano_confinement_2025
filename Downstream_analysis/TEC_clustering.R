library(readxl)
library(cluster)
library(dynamicTreeCut)


c15min = read_xlsx("NC_15m_deseq_ALL.xlsx")
colnames(c15min)[1] <- "SYMBOL"
c1h = read_xlsx("NC_1h_deseq_ALL.xlsx")
colnames(c1h)[1] <- "SYMBOL"
c4h = read_xlsx("NC_4h_deseq_ALL.xlsx")
colnames(c4h)[1] <- "SYMBOL"
c24h = read_xlsx("NC_ON_deseq_ALL.xlsx")
colnames(c24h)[1] <- "SYMBOL"


fc_data = data.frame(gene = c15min$SYMBOL, FC_c15min = c15min$log2FoldChange,
                     FC_c1h = c1h$log2FoldChange, FC_c4h = c4h$log2FoldChange,
                     fc_24h = c24h$log2FoldChange)
fc_data <- na.omit(fc_data)
rownames(fc_data) = fc_data$gene
fc_data$gene = NULL


dist_matrix = dist(fc_data)
hc <- hclust(dist_matrix, method = "ward.D2")


cluster_assignment <- cutreeDynamic(
  hc, 
  distM = as.matrix(dist_matrix), 
  deepSplit = 0,
  minClusterSize = 100,
)

gene_clusters = data.frame(Gene = rownames(fc_data), Cluster = cluster_assignment)


sil_scores <- silhouette(gene_clusters$Cluster, dist_matrix)
sil_threshold <- 0.1
noise_points <- which(sil_scores[, 3] < sil_threshold)
gene_clusters[noise_points, 2] <- 0
table(gene_clusters$Cluster)
