if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicRanges")
install.packages("remotes")
remotes::install_github("CSOgroup/CALDER2.0")

library("CALDER")
chrs = c(1:22, "X")
CALDER(contact_file_hic="file.hic", 
			chrs=chrs, bin_size=5000,genome='hg38',save_intermediate_data=FALSE,sub_domains=FALSE,
			save_dir="output")
