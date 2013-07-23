library(GenomicRanges)

source("shared_code/rdata_common.r")
source("shared_code/granges_common.r")

calc_pol_enrichments_for_txs <- function(txs) {
  piv <- txs$strand == 1
  niv <- txs$strand == -1

  txs$tss_start            <- NA
  txs$tss_end              <- NA
  txs$tu_start             <- NA
  txs$tu_end               <- NA
  txs$downstream_tss_start <- NA
  txs$downstream_tss_end   <- NA

  txs[piv, ] <- transform(txs[piv, ], tss_start      = start, 
                                      tss_end        = start + 200,
                                      
                                      tu_start       = start + 400,
                                      tu_end         = end,
                                      
                                      downstream_tss_start = start + 201,
                                      downstream_tss_end   = start + 401)

  txs[niv, ] <- transform(txs[niv, ], tss_start      = end - 200, 
                                      tss_end        = end,

                                      tu_start       = start,
                                      tu_end         = end - 400,
                                      
                                      downstream_tss_start = end - 401,
                                      downstream_tss_end   = end - 201)

  # For short genes < 600bp, TU is entire gene
  short.i <- which(txs$end - txs$start + 1 <= 600)
  txs[short.i, ] <- transform(txs[short.i, ], tu_start = start,
                                              tu_end   = end)

  base_path <- rdata_base_path()

  samples <- read.delim("samples/pol_samples.txt", stringsAsFactors=F, header=T)
  samples <- subset(samples, use == 1)

  results <- NULL

  tss.gr <- with(txs, GRanges(ranges     = IRanges(start=tss_start, end=tss_end), 
                              seqnames   = chr,
                              fb_gene_id = fb_gene_id, 
                              fb_tx_id   = fb_tx_id))

  tu.gr <- with(txs, GRanges(ranges     = IRanges(start=tu_start, end=tu_end), 
                             seqnames   = chr,
                             fb_gene_id = fb_gene_id, 
                             fb_tx_id   = fb_tx_id))

  ds.gr <- with(txs, GRanges(ranges     = IRanges(start=downstream_tss_start, end=downstream_tss_end), 
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

      values(sample.gr)$tu.ip.signal  <- regionSums(tu.gr, ip.cov)  
      values(sample.gr)$tu.wce.signal <- regionSums(tu.gr, wce.cov) 

      values(sample.gr)$downstream.tss.ip.signal  <- regionSums(ds.gr, ip.cov)
      values(sample.gr)$downstream.tss.wce.signal <- regionSums(ds.gr, wce.cov)

      sample.df <- as.data.frame(sample.gr)
      sample.df$seqnames <- NULL
      sample.df$start    <- NULL
      sample.df$end      <- NULL
      sample.df$width    <- NULL
      sample.df$strand   <- NULL

      sample.df$tp <- timepoint
      sample.df$replicate <- samples.tp$name[i]
      sample.df$antibody  <- samples.tp$antibody[i]
      sample.df$total_ip_signal <- ip.sum
      sample.df$total_wce_signal <- wce.sum
      sample.df <- transform(sample.df, tss.ratio  = log2( (tss.ip.signal/ip.sum) / (tss.wce.signal/wce.sum)),
                                        tu.ratio   = log2( (tu.ip.signal/ip.sum)  / (tu.wce.signal/wce.sum)),
                                        dst.ratio  = log2( (downstream.tss.ip.signal/ip.sum) / (downstream.tss.wce.signal/wce.sum)))

      results <- rbind(results, sample.df)
    }
  }

  pol_tss <- results

  # select the tx with the highest TSS ratio (break ties with TU ratio)

  pol_tss$selected_tx <- TRUE
  pol_tss <- pol_tss[order(pol_tss$fb_gene_id, pol_tss$tss.ratio, pol_tss$tu.ratio, decreasing=TRUE), ]
  pol_tss$selected_tx[duplicated(pol_tss$fb_gene_id)] <- FALSE
  pol_tss
}
