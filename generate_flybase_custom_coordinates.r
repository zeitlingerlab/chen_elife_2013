source("shared_code/flybase.r")

flybase.custom <- flybase_with_custom_txs()

save(flybase.custom, file="flybase/flybase.custom.RData")

