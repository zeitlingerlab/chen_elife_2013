library(plyr)

add_quantile <- function(df, signal_col) {
  new_col <- paste(signal_col, ".q", sep="")
  
  qfun <- ecdf(df[, signal_col])
  df$tmp_col <- qfun(df[, signal_col])
  names(df)[which(names(df) == "tmp_col")] <- new_col
  df
}

summarize_tss_tu_values <- function(pol_tss, all.txs) {

  genes <- transform(all.txs, tx_width = end - start + 1)
  tx.gene_map <- unique(genes[, c("fb_gene_id", "fb_symbol", "fb_tx_id", "tx_width")])

  if(length(is.na(pol_tss$tss.ratio)) > 0) pol_tss$tss.ratio[is.na(pol_tss$tss.ratio)] <- -Inf
  if(length(is.na(pol_tss$tu.ratio)) > 0)  pol_tss$tu.ratio[is.na(pol_tss$tu.ratio)] <- -Inf

  q.values <- ddply(pol_tss, .(replicate), function(x) { add_quantile(x, "tss.ip.signal")[, c("fb_tx_id", "tss.ip.signal.q")] })
  message("Merging...")
  pol_tss <- merge(pol_tss, q.values)

  pol_tss <- transform(pol_tss, tss.enriched = tss.ratio >= 1,
                                si           = tss.ratio - tu.ratio)

  tss.wide   <- reshape(pol_tss[, c("fb_tx_id", "tss.enriched", "replicate")], 
                        idvar="fb_tx_id", v.names="tss.enriched", timevar="replicate", direction="wide")
  tu.wide    <- reshape(pol_tss[, c("fb_tx_id", "tu.ratio", "replicate")],
                        idvar="fb_tx_id", v.names="tu.ratio", timevar="replicate", direction="wide")
  tss.q.wide <- reshape(pol_tss[, c("fb_tx_id", "tss.ip.signal.q", "replicate")],
                        idvar="fb_tx_id", v.names="tss.ip.signal.q", timevar="replicate", direction="wide")
  si.wide    <- reshape(pol_tss[, c("fb_tx_id", "si", "replicate")],
                        idvar="fb_tx_id", v.names="si", timevar="replicate", direction="wide")

  tss.wide$tss_count                 <- rowSums(as.matrix(tss.wide[, -1]))
  tu.wide$mean_tu                    <- log2(rowMeans(2^as.matrix(tu.wide[, -1])))
  tss.q.wide$min_tss_signal_quantile <- apply(as.matrix(tss.q.wide[, -1]), 1, min)
  si.wide$mean_si                    <- rowMeans(as.matrix(si.wide[, -1]))

  counts.df <- merge_recurse(list(tss.wide[, c("fb_tx_id", "tss_count")],
                                  tss.q.wide[, c("fb_tx_id", "min_tss_signal_quantile")],
                                  si.wide[, c("fb_tx_id", "mean_si")],
                                  tu.wide[, c("fb_tx_id", "mean_tu")]))
  merge(tx.gene_map, counts.df)
}
