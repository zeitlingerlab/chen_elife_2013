source("shared_code/flybase.r")

rna <- read.delim("rnaseq/journal.pbio.1000590.s002.txt", stringsAsFactors=F)
genes <- unique(flybase_txs()[, c("fb_gene_id", "fb_symbol")])

rna$fb_gene_id <- genes$fb_gene_id[match(rna$NAME, genes$fb_symbol)]

samples <- names(rna)[7:30]

rna.df <- NULL

for(s in samples) {
  message("Sample: ", s)
  w1  <- sprintf("w1_%s", s)
  cas <- sprintf("cas_%s", s)
  
  sample.df <- rna[, c("fb_gene_id", s, w1, cas)]
  names(sample.df)[2:4] <- c("RPKM.total", "RPKM.w1", "RPKM.cas")
  sample.df$sample <- s
  
  rna.df <- rbind(rna.df, sample.df)
}

rna.matpat <- rna.df
save(rna.matpat, file="rnaseq/rna.matpat.RData")
