# To produce barcode plot examples
# Goknur Giner
# 23 October 2016
# rm(list=ls())
library(limma); library(edgeR)

# Load data
load("DGEListFilteredNormalized.RData")

# Design matrix
Subtype <- d$samples$group
design <- model.matrix(~Subtype)
design

# voom
v <- voom(d, design, plot=TRUE)
fit <- lmFit(v,design)
fit <- eBayes(fit,robust=TRUE)
summary(decideTests(fit))

# Gene sets indices
load("human_c2_v5.rdata")
length(idx.msig <- names(Hs.c2))
# length(idx.kegg <- grep("^B", names(Hs.c2)))
msig <- Hs.c2[idx.msig]
head(msig, n = 3)
idx.msig <- ids2indices(msig, rownames(v$genes))
length(idx.msig)

# fry
system.time(fry.table <- fry(v, index = idx.msig, design = design, contrast = 2, standardize = "posterior.sd"))
head(fry.table,20)
write.table(head(fry.table,20),paste0("fry-", Sys.Date(),".txt"))
vant_up <- idx.msig[grep("VANTVEER_BREAST_CANCER_ESR1_UP", names(idx.msig))]
vant_down <- idx.msig[grep("VANTVEER_BREAST_CANCER_ESR1_DN", names(idx.msig))]
vant_up
vant_down

# roast
system.time(roast.table <- roast(v, index = idx.msig, design = design, nrot = 99999))
write.table(roast.table, paste0("roast_1e6rot-", Sys.Date(),".txt"))

# barcode plot
setEPS()
postscript(file = "..//figures//barcode_plot.eps")
barcodeplot(fit$t[,2], index = vant_up[[1]], index2 = vant_down[[1]])
dev.off()





