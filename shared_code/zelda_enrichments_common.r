library(GenomicRanges)
library(plyr)

source("shared_code/granges_common.r")
source("shared_code/rdata_common.r")

calc_zelda_enrichments_for_txs <- function(txs) {
  piv <- txs$strand == 1
  niv <- txs$strand == -1

  txs$promoter_start      <- NA
  txs$promoter_end        <- NA

  txs[piv, ] <- transform(txs[piv, ], promoter_start      = start - 2000, 
                                      promoter_end        = start)

  txs[niv, ] <- transform(txs[niv, ], promoter_start      = end, 
                                      promoter_end        = end + 2000)

  base_path <- "~/rdata/chromatin_early_embryo/"

  tss.gr <- with(txs, GRanges(ranges     = IRanges(start=promoter_start, end=promoter_end), 
                              seqnames   = chr,
                              fb_gene_id = fb_gene_id, 
                              fb_tx_id   = fb_tx_id))
  results <- NULL

  # Zelda nc8 data from Eisen lab
  ip.cov  <- get_load("eisen_zelda_chipseq/zld_cycle_8.cov.RData")
  wce.cov <- get_load(paste0(rdata_base_path(), "preMBT_wce_1.cov.RData"))

  ip.sum  <- total_signal(ip.cov)
  wce.sum <- total_signal(wce.cov)

  sample.gr <- tss.gr

  values(sample.gr)$promoter.ip.signal  <- regionSums(tss.gr, ip.cov)
  values(sample.gr)$promoter.wce.signal <- regionSums(tss.gr, wce.cov)

  sample.df <- as.data.frame(sample.gr)
  sample.df$seqnames <- NULL
  sample.df$start    <- NULL
  sample.df$end      <- NULL
  sample.df$width    <- NULL
  sample.df$strand   <- NULL

  sample.df$total_ip_signal <- ip.sum
  sample.df$total_wce_signal <- wce.sum
  sample.df <- transform(sample.df, zelda.ratio = log2( (promoter.ip.signal/ip.sum) / (promoter.wce.signal/wce.sum)))
  results <- rbind(results, sample.df)
  results
}
