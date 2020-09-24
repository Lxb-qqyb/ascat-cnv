args<-commandArgs(T)
path<-args[1]
setwd(path)
library(ASCAT)
ascat.bc = ascat.loadData("Tumor_LogR.txt","Tumor_BAF.txt","Germline_LogR.txt","Germline_BAF.txt")
ascat.plotRawData(ascat.bc)
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)
saveRDS(ascat.output,"ascat.output.rds")
segTab<-ascat.output$segments
write.table(file="segment.txt",segTab)
ploidy<-ascat.output$ploidy
write.table(file="ploidy.txt",ploidy)
