source("shared_code/zelda_enrichments_common.r")
source("shared_code/flybase.r")

txs <- flybase_with_custom_txs()

zelda_tss <- calc_zelda_enrichments_for_txs(txs)

save(zelda_tss, file="zelda_tss.custom.RData")



