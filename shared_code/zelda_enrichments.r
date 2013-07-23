library(ggplot2)
library(scales)
library(parallel)

source("shared_code/flybase.r")
source("shared_code/load_groups.r")

pn <- function(i) {
  prettyNum(i, big.mark=",")
}

zelda <- subset(get(load("zelda_motif_counts.RData")), mm == 0)
fb.txs <- flybase_with_custom_txs()

zelda_enrichments_random_sampling <- function(zelda.counts, all.txs, group.size, count) {
  all.txs <- merge(all.txs, zelda.counts)
  all.txs <- transform(all.txs, tss = ifelse(strand == 1, start, end))

  uniq.txs <- unique(all.txs[, c("chr", "tss", "strand", "count")])
  
  total.zelda <- sum(uniq.txs$count)
  exp <- total.zelda/nrow(uniq.txs)
  
  results <- c()
  
  sample_and_calc_zelda_enrichment <- function(all.txs, group.size, expected) {
    group.uniq.txs <- unique(all.txs[sample(nrow(all.txs), group.size), c("chr", "tss", "strand", "count")])
    obs <- as.numeric(sum(group.uniq.txs$count)/nrow(group.uniq.txs))
    as.numeric(obs/expected)
  }
  
  results <- unlist(mclapply(1:count, function(x) sample_and_calc_zelda_enrichment(all.txs, group.size, exp), mc.cores=4))
  results
}

zelda_enrichments <- function(zelda.counts, all.txs, tx.ids, group_name) {
  message("Zelda enrichment for: ", group_name)
  
  all.txs <- merge(all.txs, zelda.counts)
  all.txs <- transform(all.txs, tss = ifelse(strand == 1, start, end))

  uniq.txs <- unique(all.txs[, c("chr", "tss", "strand", "count")])
  
  total.zelda <- sum(uniq.txs$count)
  
  message(" ", pn(nrow(uniq.txs)), " unique promoters")
  message(" ", pn(total.zelda), " Zelda motifs")
  message(" Expected: ", round(total.zelda/nrow(uniq.txs), 2), " motifs per promoter")
  
  group.txs <- subset(all.txs, fb_tx_id %in% tx.ids)
  group.uniq.txs <- unique(group.txs[, c("chr", "tss", "strand", "count")])
  
  group.total <- sum(group.uniq.txs$count)
  
  message(" Group has ", pn(nrow(group.uniq.txs)), " unique promoters")
  message(" Group has ", pn(group.total), " Zelda motifs")
  message(" ", round(group.total/nrow(group.uniq.txs), 2), " motifs per promoter")
  
  obs <- group.total/nrow(group.uniq.txs)
  exp <- total.zelda/nrow(uniq.txs)
  
  enrichment <- obs / exp
  message(" - enrichment: ", round(enrichment, 2))
  
  message(" - randomly sampling for ", length(tx.ids), " transcripts 10,000 times...")
  random.10k <- zelda_enrichments_random_sampling(zelda.counts, all.txs, length(tx.ids), 10000)
  
  pv <- length(which(random.10k > enrichment)) / length(random.10k)
  
  message(" - pvalue: ", round(pv, 5))
  
  data.frame(enrichment=enrichment, pvalue=pv, group_name=group_name)
}


