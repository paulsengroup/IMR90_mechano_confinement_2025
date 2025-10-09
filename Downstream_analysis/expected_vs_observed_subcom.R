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

#Plot
pheatmap(agg_tbl_matri, main = "", cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 0
         ,fontsize = 12, display_numbers = round(agg_tbl_matri, digits = 3)) 

setHook("grid.newpage", NULL, "replace")
grid.text("NC", y=-0.02, gp=gpar(fontsize=16))
grid.text("C", x=-0.02, rot=90, gp=gpar(fontsize=16))


