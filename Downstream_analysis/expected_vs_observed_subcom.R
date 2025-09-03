library(pheatmap)
library(grid)


#Subcompartment NC vs C count
sub_com = read.csv("NC_vs_C_cons_2.tsv", sep = "\t")
colnames(sub_com) = c("A", "B", "count")

#Subcompartment NC vs C in Mbp
sub_com$count = sub_com$count*5000/1000000
#Subcompartment NC vs C convert to table
agg_tbl_matri = xtabs(count ~ A + B, sub_com)
agg_tbl_matri = agg_tbl_matri[c(8, 7, 6, 5, 1:4), c(8, 7, 6, 5, 1:4)]


#subcom share in C
a = rowSums(agg_tbl_matri)/sum(rowSums(agg_tbl_matri))
#subcom share in NC
b = colSums(agg_tbl_matri)/sum(colSums(agg_tbl_matri))

#Calculating the expected matrix based on subcompartment share in each condition
matrix = c()
for(value1 in a){
  row = c()
  for(value2 in b){
    row = append(row, value1*value2)
  }
  matrix = append(matrix, row)
}
matrix = matrix(matrix, nrow = 8, ncol=8)
rownames(matrix) = rownames(agg_tbl_matri) 
colnames(matrix) = colnames(agg_tbl_matri)
#Calculating the share of each subcompartment NC vs. C from the observed matrix 
obs = agg_tbl_matri/sum(agg_tbl_matri)
#log2 subcompartment matrix NC vs C, observed/expected matrix
ods = log2(obs/matrix)
