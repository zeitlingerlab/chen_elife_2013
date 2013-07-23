library(GenomicRanges)
library(plyr)

source("shared_code/rdata_common.r")
source("shared_code/granges_common.r")

calc_tbp_enrichments_for_txs <- function(txs) {
  piv <- txs$strand == 1
  niv <- txs$strand == -1

  txs$tss_start      <- NA
  txs$tss_end        <- NA
  txs$tu_start       <- NA
  txs$tu_end         <- NA

  txs[piv, ] <- transform(txs[piv, ], tss_start      = start - 100, 
                                      tss_end        = start)

  txs[niv, ] <- transform(txs[niv, ], tss_start      = end, 
                                      tss_end        = end + 100)

  base_path <- rdata_base_path()

  samples <- read.delim("samples/tbp_samples.txt", stringsAsFactors=F, header=T)
  samples <- subset(samples, use == 1)

  results <- NULL

  tss.gr <- with(txs, GRanges(ranges     = IRanges(start=tss_start, end=tss_end), 
                              seqnames   = chr,
                              fb_gene_id = fb_gene_id, 
                              fb_tx_id   = fb_tx_id))
  results <- NULL

  for(timepoint in sort(unique(samples$tp))) {
    message("Time: ", timepoint)

    samples.tp <- subset(samples, tp == timepoint)

    for(i in 1:nrow(samples.tp)) {
      nothing <- gc()
      message(" - sample: ", samples.tp$sample[i])
      ip_cov_file  <- paste(base_path, samples.tp$sample[i],  ".cov.RData", sep="")
      wce_cov_file <- paste(base_path, samples.tp$control[i], ".cov.RData", sep="")

      ip.cov  <- get_load(ip_cov_file)
      wce.cov <- get_load(wce_cov_file)

      ip.sum  <- total_signal(ip.cov)
      wce.sum <- total_signal(wce.cov)

      sample.gr <- tss.gr

      values(sample.gr)$tss.ip.signal  <- regionSums(tss.gr, ip.cov)
      values(sample.gr)$tss.wce.signal <- regionSums(tss.gr, wce.cov)

      sample.df <- as.data.frame(sample.gr)
      sample.df$seqnames <- NULL
      sample.df$start    <- NULL
      sample.df$end      <- NULL
      sample.df$width    <- NULL
      sample.df$strand   <- NULL

      sample.df$tp <- timepoint
      sample.df$replicate <- samples.tp$name[i]
      sample.df$total_ip_signal <- ip.sum
      sample.df$total_wce_signal <- wce.sum
      sample.df <- transform(sample.df, tss.ratio = log2( (tss.ip.signal/ip.sum) / (tss.wce.signal/wce.sum)))
      results <- rbind(results, sample.df)
    }
  }

  tbp_tss <- results

  tbp_tss
}
