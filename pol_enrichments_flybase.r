source("shared_code/pol_enrichments_common.r")
source("shared_code/flybase.r")

txs <- flybase_txs()

pol_tss <- calc_pol_enrichments_for_txs(txs)

save(pol_tss, file="pol_tss.flybase.RData")



