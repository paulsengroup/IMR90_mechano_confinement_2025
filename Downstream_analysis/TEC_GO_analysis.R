library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(cluster)
library(GOSemSim)
library(igraph)

#GO analysis

rna_clusters = read.table("clusters.tsv", header  = TRUE)

entrezids <- AnnotationDbi::select(org.Hs.eg.db, 
                                   keys = rna_clusters$Gene,
                                   columns = c("ENTREZID"),
                                   keytype = "SYMBOL")
entrezids <- entrezids[!is.na(entrezids$ENTREZID), ]
colnames(entrezids)[1] = "Gene"
rna_clusters <- merge(entrezids, rna_clusters, by = "Gene")

sample <- list(rna_clusters[rna_clusters$Cluster == "0", 2],
               rna_clusters[rna_clusters$Cluster == "1", 2],
               rna_clusters[rna_clusters$Cluster == "2", 2],
               rna_clusters[rna_clusters$Cluster == "3", 2],
               rna_clusters[rna_clusters$Cluster == "4", 2],
               rna_clusters[rna_clusters$Cluster == "5", 2],
               rna_clusters[rna_clusters$Cluster == "6", 2],
               rna_clusters[rna_clusters$Cluster == "7", 2])


names(sample) <- c("Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7")


compGO <- compareCluster(geneCluster = sample, 
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db,
                         pvalueCutoff  = 0.05, 
                         pAdjustMethod = "BH", ont = "BP")

for(i in 1:length(compGO@compareClusterResult$geneID)){
  a <- compGO@compareClusterResult$geneID[i]
  a <- strsplit(a, "/")
  a <- a[[1]]
  a <- data.frame(a)
  colnames(a) <- "ENTREZID"
  a <- merge(a, entrezids, by="ENTREZID")
  a <- paste(a$Gene, collapse='/')
  compGO@compareClusterResult$geneID[i] <- a
}


#GO similarity calculation
d <- godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)


df = NULL
for(i in 1:7){
  all_cluster_genes = read_xlsx("results.xlsx")
  all_cluster_genes <- all_cluster_genes %>%
    filter(as.numeric(sub("/.*", "", GeneRatio)) >= 10)
  
  
  cluster_name = paste("Cluster", i)
  number = 50
  all_cluster_genes = all_cluster_genes[all_cluster_genes$Cluster == cluster_name, ]
  n <- min(length(rownames(all_cluster_genes[, 2])), number)
  all_cluster_genes = all_cluster_genes[1:n, ]
  
  
  similarity_matrix <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      similarity_matrix[i, j] <- goSim(all_cluster_genes[, 2]$ID[i], 
                                       all_cluster_genes[, 2]$ID[j],
                                       semData=d)
      similarity_matrix[j, i] <- similarity_matrix[i, j]
    }
  }
  
  threshold <- 0.3
  adj_matrix <- similarity_matrix > threshold 
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
  
  components <- components(g)
  
  GO_Term <- split(all_cluster_genes[, 3], components$membership)
  GO_Term <- unlist(lapply(GO_Term, function(group) head(group, 1)), use.names = FALSE)
  
  padj <- split(all_cluster_genes[, 7], components$membership)
  padj <- unlist(lapply(padj, function(group) head(group, 1)), use.names = FALSE)
  
  GeneRatio <- split(all_cluster_genes[, 4], components$membership)
  GeneRatio <- unlist(lapply(GeneRatio, function(group) head(group, 1)), use.names = FALSE)
  Cluster = rep(cluster_name, length(padj))
  
  if(i == 1){
    df <- data.frame(
      Cluster = Cluster,
      GO_Term = GO_Term, 
      padj = padj,
      GeneRatio = GeneRatio
    )
  }else{
    df_tmp <- data.frame(
      Cluster = Cluster,
      GO_Term = GO_Term, 
      padj = padj,
      GeneRatio = GeneRatio
    )
    
    df <- rbind(df, df_tmp)
  }
  
}

df
