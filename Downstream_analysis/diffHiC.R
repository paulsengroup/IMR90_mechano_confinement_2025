library(diffHic)
library(edgeR)
library(dplyr)



group <- factor(c(
  rep("NC",3),
  rep("C15m",3),
  rep("C4h",3)
))

keep <- aveLogCPM(asDGEList(data)) > 0
data <- data[keep,]


data <- normOffsets(data, se.out=TRUE, method = "loess")
y <- asDGEList(data)

# Defining contrast groups NC vs. C15m and NC vs. C4h
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

y <- estimateDisp(y, design)

contrast_C15m_vs_NC <- makeContrasts(C15m - NC, levels = design)
contrast_C4h_vs_NC <- makeContrasts(C4h - NC, levels = design)
                   
#fit <- glmFit(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)

res_C15m <- glmQLFTest(fit, contrast = contrast_C15m_vs_NC)
res_C4h  <- glmQLFTest(fit, contrast = contrast_C4h_vs_NC)

#res_C15m <- glmLRT(fit, contrast = contrast_C15m_vs_NC)
#res_C4h  <- glmLRT(fit, contrast = contrast_C4h_vs_NC)


tab_C15m <- res_C15m$table
tab_C4h  <- res_C4h$table

tab_C15m$FDR <- p.adjust(tab_C15m$PValue, method="BH")
tab_C4h$FDR  <- p.adjust(tab_C4h$PValue, method="BH")

