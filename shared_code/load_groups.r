suppressPackageStartupMessages(library(data.table))

# First wave genes (pre-MBT)

load_fw_groups <- function() {
  fw.df <- get(load("step4_final_first_wave_gene_list_output/fw.df.RData"))
  cr2012 <- read.delim("gaertner_johnston_2012/polii_tss_enrichments_mef2_6to8h.txt", stringsAsFactors=FALSE, header=TRUE)
  
  # Paused first wave genes: those with mean TU enrichment below 0 (log2)
  fw.paused <- subset(fw.df, mean_tu < 0)
  
  fw.not_paused <- subset(fw.df, mean_tu >= 0)

  # Paused later: genes with TSS enrichment in the top 20% in 6-8h muscle time course
  rep1.top <- quantile(cr2012$TSS.rep1.mef2_06to08h, 0.80, na.rm=TRUE)
  rep2.top <- quantile(cr2012$TSS.rep2.mef2_06to08h, 0.80, na.rm=TRUE)
  cr2012 <- transform(cr2012, TSS.enriched = TSS.rep1.mef2_06to08h > rep1.top & TSS.rep2.mef2_06to08h > rep2.top)
  high_tss_genes.6 <- as.character(subset(cr2012, TSS.enriched == TRUE)$flybase_gene_id)

  fw.paused_later <- subset(fw.not_paused, fb_gene_id %in% high_tss_genes.6)
  
  fw.never_paused <- subset(fw.not_paused, ! fb_gene_id %in% high_tss_genes.6)

  list("paused"=fw.paused, 
       "paused_later"=fw.paused_later,
       "never_paused"=fw.never_paused)
}

fw_groups <- load_fw_groups()
fw_all    <- do.call(rbind, fw_groups)

# Second wave genes (MBT)

load_sw_groups <- function() {
  sw.df <- get(load("step6_select_second_wave_genes_by_tss_output/sw.df.RData"))

  rna.df <- get(load("rnaseq/eisen.rnaseq.RData"))
  
  maternal.exp <- as.character(subset(rna.df, sample == "cc10" & RPKM >= 1)$fb_gene_id)
  high.exp     <- as.character(subset(rna.df, sample == "cc14d" & RPKM >= 5)$fb_gene_id)
  
  sw.maternal.df <- subset(sw.df, fb_gene_id %in% maternal.exp)
  sw.dev.df      <- subset(sw.df, ! fb_gene_id %in% maternal.exp)
  sw.dev.high.df <- subset(sw.dev.df, fb_gene_id %in% high.exp)
  sw.dev.low.df  <- subset(sw.dev.df, ! fb_gene_id %in% high.exp)

  list("maternal"=sw.maternal.df, 
       "dev_high"=sw.dev.high.df,
       "dev_low"=sw.dev.low.df)
}

sw_groups <- load_sw_groups()
sw_all <- do.call(rbind, sw_groups)

# Inactive genes

load_inactive_genes <- function() {
  insitu <- unique(as.data.table(get(load("imago/insitu.RData"))[, c("fb_gene_id", "go_term", "stage_name")]))
  no_staining.dt <- insitu[, list(no_staining_count=length(which(go_term == "no staining"))), by=fb_gene_id]
  no_staining_ids <- unique(as.character(subset(no_staining.dt, no_staining_count == 6)$fb_gene_id))
  
  no_staining_ids[!no_staining_ids %in% c(fw_all$fb_gene_id, sw_all$fb_gene_id)]
}

inactive_gene_ids <- load_inactive_genes()

