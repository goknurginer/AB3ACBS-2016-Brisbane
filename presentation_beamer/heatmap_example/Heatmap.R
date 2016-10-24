# Script for Heatmap for the useR 2015 conference presentation
# Author: Goknur Giner 
# Date created:19 June 2016
# Date modified: 19 June 2016
setwd("C:/Users/giner.g.WEHI/Dropbox/useR2015_presentation/Presentation/Heatmap")
rm(list = ls());library(pheatmap)

library(pheatmap)
# Generate some data
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = LETTERS[1:10]
rownames(test) = paste0("Gene", 1:20)
# original figure
pheatmap(test)
annotation <- data.frame(Condition = factor(1:10 %% 2 == 0, labels = c("Experiment1", "Experiment2")))
rownames(annotation) <- colnames(test) # check out the row names of annotation

pheatmap(test, annotation = annotation, filename = "Heatmap.pdf")


setEPS()
postscript("Heatpmap.eps")
pheatmap(test, annotation = annotation)
dev.off()
