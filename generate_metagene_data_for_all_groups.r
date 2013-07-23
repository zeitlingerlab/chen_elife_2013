library(parallel)
library(ggplot2)
library(GenomicRanges)

source("shared_code/metagene_graphs.R")
source("shared_code/flybase.r")
source("shared_code/load_groups.r")
source("shared_code/rdata_common.r")

tp2.pol.r1 <- "preMBT_pol_1"
tp2.pol.r2 <- "preMBT_pol_2"
tp2.pol.r3 <- "preMBT_pol_3"
tp2.pol.r4 <- "preMBT_pol_4"

tp2.wce.r1 <- "preMBT_wce_1"

tp2.tbp.r1 <- "preMBT_tbp_1"
tp2.tbp.r2 <- "preMBT_tbp_2"

tp2.k4.r1 <- "preMBT_k4_1"
tp2.k4.r2 <- "preMBT_k4_2"

tp3.pol.r1 <- "MBT_pol_1"
tp3.pol.r2 <- "MBT_pol_2"
tp3.pol.r3 <- "MBT_pol_3"

tp3.tbp.r1 <- "MBT_tbp_1"
tp3.tbp.r2 <- "MBT_tbp_2"

tp3.k4.r1 <- "MBT_k4_1"
tp3.k4.r2 <- "MBT_k4_2"

tp3.h3ac.r1 <- "MBT_h3ac_1"
tp3.h3ac.r2 <- "MBT_h3ac_2"

tp3.wce.r1 <- "MBT_wce_1"

samples.df <- NULL

add_sample <- function(samples.df, ip, wce, tp, name, replicate) {
  rbind(samples.df, data.frame(stringsAsFactors=F, name=name, replicate=replicate, tp=tp, ip=ip, wce=wce))
}

samples.df <- add_sample(samples.df, tp2.pol.r1, tp2.wce.r1, "pre-MBT", "Pol II", 1)
samples.df <- add_sample(samples.df, tp2.pol.r2, tp2.wce.r1, "pre-MBT", "Pol II", 2)
samples.df <- add_sample(samples.df, tp2.pol.r3, tp2.wce.r1, "pre-MBT", "Pol II", 3)
samples.df <- add_sample(samples.df, tp2.pol.r4, tp2.wce.r1, "pre-MBT", "Pol II", 4)

samples.df <- add_sample(samples.df, tp2.tbp.r1, tp2.wce.r1, "pre-MBT", "TBP", 1)
samples.df <- add_sample(samples.df, tp2.tbp.r2, tp2.wce.r1, "pre-MBT", "TBP", 2)

samples.df <- add_sample(samples.df, tp2.k4.r1, tp2.wce.r1, "pre-MBT", "K4", 1)
samples.df <- add_sample(samples.df, tp2.k4.r2, tp2.wce.r1, "pre-MBT", "K4", 2)

samples.df <- add_sample(samples.df, tp3.pol.r1, tp3.wce.r1, "MBT", "Pol II", 1)
samples.df <- add_sample(samples.df, tp3.pol.r2, tp3.wce.r1, "MBT", "Pol II", 2)
samples.df <- add_sample(samples.df, tp3.pol.r3, tp3.wce.r1, "MBT", "Pol II", 3)

samples.df <- add_sample(samples.df, tp3.tbp.r1, tp3.wce.r1, "MBT", "TBP", 1)
samples.df <- add_sample(samples.df, tp3.tbp.r2, tp3.wce.r1, "MBT", "TBP", 2)

samples.df <- add_sample(samples.df, tp3.k4.r1, tp3.wce.r1, "MBT", "K4", 1)
samples.df <- add_sample(samples.df, tp3.k4.r2, tp3.wce.r1, "MBT", "K4", 2)

samples.df <- add_sample(samples.df, tp3.h3ac.r1, tp3.wce.r1, "MBT", "H3Ac", 1)
samples.df <- add_sample(samples.df, tp3.h3ac.r2, tp3.wce.r1, "MBT", "H3Ac", 2)

all.txs <- flybase_with_custom_txs()
names(all.txs)[which(names(all.txs) == "end")] <- "stop"

gene_groups <- list("FW All"=fw_all$fb_tx_id, 
                    "FW Paused"=fw_groups$paused$fb_tx_id,
                    "FW Paused later"=fw_groups$paused_later$fb_tx_id,
                    "FW Not paused"=fw_groups$never_paused$fb_tx_id,
                    "SW Maternal"=sw_groups$maternal$fb_tx_id,
                    "SW Dev high"=sw_groups$dev_high$fb_tx_id,
                    "SW Dev low"=sw_groups$dev_low$fb_tx_id)

cov_cache <- list()
load_and_cache_cov <- function(obj_name) {
  if(is.null(cov_cache[[obj_name]])) {
    cov_cache[[obj_name]] <<- get_load(paste0(rdata_base_path(), "/", obj_name, ".cov.RData"))
  }
  cov_cache[[obj_name]]
}

# Preload cache

for(i in 1:nrow(samples.df)) {
  ip.cov <- load_and_cache_cov(samples.df$ip[i])
  wce.cov <- load_and_cache_cov(samples.df$wce[i])
}
nothing <- gc()
rm(ip.cov)
rm(wce.cov)

parallel_metagene_for_gene_group <- function(i, samples.df, group.txs) {
  message(" - ", samples.df$name[i], " tp ", samples.df$tp[i], " replicate ", samples.df$replicate[i])
  reads.loop <- get_metagene_enrichment(load_and_cache_cov(samples.df$ip[i]), 
                                        load_and_cache_cov(samples.df$wce[i]),
                                        group.txs, before_tss=250, after_tss=1000, smooth=25, sample=samples.df$name[i])
  reads.loop$replicate <- samples.df$replicate[i]
  reads.loop$tp <- samples.df$tp[i]
  reads.loop
}

metagene.df <- NULL
for(g in names(gene_groups)) {
  nothing <- gc()
  message("Gene group: ", g)
  group.txs <- subset(all.txs, fb_tx_id %in% gene_groups[[g]])
  reads.group <- mclapply(1:nrow(samples.df), function(i) { parallel_metagene_for_gene_group(i, samples.df, group.txs) }, mc.cores=2)
  reads.group <- do.call(rbind, reads.group)
  reads.group$gene_group <- g
  reads.group$group_size <- nrow(group.txs)
  metagene.df <- rbind(metagene.df, reads.group)
  }

save(metagene.df, file="metagene.df.RData")
