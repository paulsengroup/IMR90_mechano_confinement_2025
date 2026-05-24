library(multiHiCcompare)


hicexp_c15m <- make_hicexp(c15min_rep1, c15min_rep2, c15min_rep3, nc_rep1, nc_rep2, nc_rep3,
                           groups = c("C15m", "C15m","C15m", "NC", "NC","NC"))

hicexp_c15m <- cyclic_loess(hicexp_c15m)

hicexp_c15m <- hic_exactTest(hicexp_c15m)

res_c15m = hicexp_c15m@comparison
res_sig_c15m = res_c15m[res_c15m$p.adj < 0.05, ]


hicexp_c4h <- make_hicexp(c4h_rep1, c4h_rep2, c4h_rep3, nc_rep1, nc_rep2, nc_rep3,
                          groups = c("C4h", "C4h","C4h", "NC", "NC","NC"))

hicexp_c4h <- cyclic_loess(hicexp_c4h)

hicexp_c4h <- hic_exactTest(hicexp_c4h)

res_c4h = hicexp_c4h@comparison
res_sig_c4h = res_c4h[res_c4h$p.adj < 0.05, ]
