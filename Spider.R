## myscript.R

library(spider)
myArgs <- commandArgs(trailingOnly=TRUE)
fasta<-read.dna(myArgs, format="fasta")
Dist <- dist.dna(fasta, pairwise.deletion = TRUE)
Thresh <- localMinima(Dist)
OT <- Thresh$localMinima[1] *1
OT2 <- Thresh$localMinima[1] *100
OT2
spp <- dimnames(fasta)[[1]]
Clust <- tclust(Dist, threshold=OT)
list <- lapply(Clust, function(x) spp[x])
list