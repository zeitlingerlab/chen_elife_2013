library(reshape)
library(ggplot2)

cuff.4 <- read.delim("we04to06h_cufflinks/genes.fpkm_tracking.gz", stringsAsFactors=F, header=T)
cuff.6 <- read.delim("we06to08h_cufflinks/genes.fpkm_tracking.gz", stringsAsFactors=F, header=T)

reformat_cufflinks <- function(df, label) {
  df <- subset(df, FPKM_status == "OK")[, c("gene_id", "FPKM")]
  names(df) <- c("fb_gene_id", "RPKM")
  df$sample <- label
  df
}

we.rnaseq <- rbind(reformat_cufflinks(cuff.4, "we04to06h"), reformat_cufflinks(cuff.6, "we06to08h"))

save(we.rnaseq, file="we.rnaseq.RData")
