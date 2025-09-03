"""
#subsampling
cooltools random-sample C4h_merged.mcool::/resolutions/1000  C4h_merged_subsample_1kbp.cool -c 1205642737
cooltools random-sample NC_merged.mcool::/resolutions/1000  NC_merged_subsample_1kbp.cool -c 1205642737
cooler balance  NC_merged_subsample_1kbp.cool
cooler balance C4h_merged_subsample_1kbp.cool

cooler zoomify -r 1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000 C4h_merged_subsample_1kbp.cool -o C4h_merged_subsample.mcool
cooler zoomify -r 1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000 NC_merged_subsample_1kbp.cool -o NC_merged_subsample.mcool
"""

import numpy as np
import pandas as pd
import csv
import hictkpy

c = hictkpy.File('C_sub_sample.mcool', resolution = 50000)

chroms = list(c.chromosomes().keys())
all_counts = []
for chr_first in chroms:
    intermed_counts = []
    for chr_second in chroms:
        count = c.fetch(chr_first, chr_second, normalization="weight")
        intermed_counts.append(np.nanmean(count.to_numpy()))
    all_counts.append(intermed_counts)


with open("C_sub_sample.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(all_counts)




