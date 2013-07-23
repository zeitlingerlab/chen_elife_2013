source("shared_code/tbp_enrichments_common.r")
source("shared_code/flybase.r")

txs <- flybase_with_custom_txs()

tbp_tss <- calc_tbp_enrichments_for_txs(txs)

save(tbp_tss, file="tbp_tss.custom.RData")



