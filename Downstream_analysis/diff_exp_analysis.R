#Author: Oda Hovet


library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DESeq2)
library('org.Hs.eg.db')
library(writexl)
library(dplyr)


# Load counts table
counts <- read.table(file = "salmon_gene_counts/salmon.merged.gene_counts.tsv", sep = '\t', header = TRUE)

# Filter genes using hg38 annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_ids <- genes(txdb)@elementMetadata@listData[["gene_id"]]

# Map ENTREZ IDs to gene symbols
gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys = gene_ids,
                                      columns = "SYMBOL",
                                      keytype = "ENTREZID") %>%
  filter(!is.na(SYMBOL))


# Keep only genes present in hg38
countData <- counts %>%
  filter(gene_name %in% gene_symbols$SYMBOL) %>%
  select(-gene_id) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), sum)) %>%
  as.data.frame()


# Set gene names as row names
rownames(countData) <- countData$gene_name
countData$gene_name <- NULL
countData <- countData[, !(colnames(countData) %in% c("gene_id", "15m_REP2"))]


# Define batch variable
batch_ex_deseq <- c(2, 2, 2, 2, 2, 1, 
                    2, 2, 2, 2, 2, 2, 
                    1, 1, 1, 1, 1, 1, 
                    1, 1, 1, 1, 1, 1, 
                    2, 2, 2, 2, 1, 1, 1)
batch_ex_deseq <- factor(batch_ex_deseq)

# Column metadata 
col_data_ex_deseq <- data.frame(
  Sample = colnames(countData),
  Condition = c(rep("15m", 6),
                rep("X1h", 6),
                rep("X4h", 6),
                rep("NC", 7),
                rep("ON", 6))
)


# Add batch variable
col_data_ex_deseq$batch <- batch_ex_deseq

# DESeq2 comparison function
run_deseq_comparison <- function(count_data, col_data, 
                                 groupA, groupB, outfile_prefix) {
  # Subset data for the two conditions
  subset_idx <- which(col_data$Condition %in% c(groupA, groupB))
  counts_sub <- count_data[, subset_idx]
  meta_sub <- col_data[subset_idx, ]
  
  # Build DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_sub),
    colData = meta_sub,
    design = ~ batch + Condition
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Filter significant results
  res_filtered <- res[!is.na(res$padj) & res$padj < 0.05, ]
  res_filtered <- as.data.frame(res_filtered)
  
  # Write results to Excel
  out_file <- paste0("corrected_resfiltered_", groupA, "_vs_", groupB, "_ex_deseq.xlsx")
  write.xlsx(res_filtered, out_file, rowNames = TRUE)
  
  return(res_filtered)
}

comparisons <- list(
  c("NC", "15m"),
  c("NC", "X1h"),
  c("NC", "X4h"),
  c("NC", "ON")
)

results_list <- lapply(comparisons, function(x) {
  run_deseq_comparison(countData, col_data_ex_deseq, x[1], x[2], "ex_deseq")
})

names(results_list) <- sapply(comparisons, function(x) paste(x, collapse = "_vs_"))



results_list$NC_vs_15m
results_list$NC_vs_X1h
results_list$NC_vs_X4h
results_list$NC_vs_ON




