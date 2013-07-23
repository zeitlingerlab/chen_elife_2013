library(GenomicRanges)

flybase_txs <- function() {
  txs <- get(load("flybase/fb.transcripts.r5.47.RData"))
  valid.chrs <- c("chr2L", "chr2LHet", "chr2R", "chr2RHet", "chr3L", "chr3LHet", "chr3R", "chr3RHet", "chr4", "chrX", "chrXHet", "chrYHet")
  txs <- subset(txs, chr %in% valid.chrs)
  txs
}

flybase_txs_granges <- function() {
  genes <- flybase_txs()
  genes.gr <- with(genes, GRanges(ranges     = IRanges(start=start, end=end), 
                                  seqnames   = chr,
                                  strand     = ifelse(strand == 1, "+", "-"), 
                                  fb_tx_id   = fb_tx_id,
                                  fb_gene_id = fb_gene_id, 
                                  fb_symbol  = fb_symbol))
  genes.gr
}

flybase_with_custom_txs <- function() {
  custom.txs <- get(load("step2_define_custom_transcripts_output/custom_txs.RData"))
  rbind(custom.txs, flybase_txs())
}

flybase_with_custom_txs_granges <- function() {
  genes <- flybase_with_custom_txs()
  genes.gr <- with(genes, GRanges(ranges     = IRanges(start=start, end=end), 
                                  seqnames   = chr,
                                  strand     = ifelse(strand == 1, "+", "-"), 
                                  fb_tx_id   = fb_tx_id,
                                  fb_gene_id = fb_gene_id, 
                                  fb_symbol  = fb_symbol))
  genes.gr
}

