library(reshape)
library(ggplot2)

cuffdiff <- read.delim("eisen_cuffdiff/genes.fpkm_tracking.gz", stringsAsFactors=F, header=T)

fpkm.columns <- grep("_FPKM$", names(cuffdiff))

cuffdiff.m <- melt(cuffdiff[, c(1, fpkm.columns)])
names(cuffdiff.m) <- c("fb_gene_id", "sample", "RPKM")
cuffdiff.m$sample <- gsub("(.*)_FPKM$", "\\1", cuffdiff.m$sample)

eisen.rnaseq <- cuffdiff.m

save(eisen.rnaseq, file="eisen.rnaseq.RData")

library(ggplot2)

g <- ggplot(eisen.rnaseq, aes(x=log2(RPKM+1), color=sample)) +
     geom_line(stat="density") + 
     theme_bw() 

