library(GenomicRanges)
library(lattice)

get_load <- function(filename) {
  message("Loading: ", filename)
  get(load(filename))
}

source("shared_code/flybase.r")
source("shared_code/heatmap.r")
source("shared_code/load_groups.r")
source("shared_code/rdata_common.r")

genes <- flybase_with_custom_txs()

genes.fw        <- subset(genes, fb_tx_id %in% fw_all$fb_tx_id)
genes.sw.mat    <- subset(genes, fb_tx_id %in% sw_groups$maternal$fb_tx_id)
genes.sw.dev.lo <- subset(genes, fb_tx_id %in% sw_groups$dev_low$fb_tx_id)
genes.sw.dev.hi <- subset(genes, fb_tx_id %in% sw_groups$dev_high$fb_tx_id)

message("1st wave count: ", nrow(genes.fw))
message("2nd wave maternal: ", nrow(genes.sw.mat))
message("2nd wave dev low count:  ", nrow(genes.sw.dev.lo))
message("2nd wave dev high count: ", nrow(genes.sw.dev.hi))

cache_file <- "figure_1c_heatmap_reads.cache.RData"
if(file.exists(cache_file)) {
  message("Loading enrichment reads from cache...")
  load(cache_file)
} else {
  message("Generating enrichment reads...")
 
  pol.2h <- get_load(paste0(rdata_base_path(), "enrichment/preMBT_pol_minimum_enrichment.cov.RData"))
  tbp.2h <- get_load(paste0(rdata_base_path(), "enrichment/preMBT_tbp_minimum_enrichment.cov.RData"))
  k4.2h  <- get_load(paste0(rdata_base_path(), "enrichment/preMBT_k4_minimum_enrichment.cov.RData"))

  pol.3h  <- get_load(paste0(rdata_base_path(), "enrichment/MBT_pol_minimum_enrichment.cov.RData"))
  tbp.3h  <- get_load(paste0(rdata_base_path(), "enrichment/MBT_tbp_minimum_enrichment.cov.RData"))
  k4.3h   <- get_load(paste0(rdata_base_path(), "enrichment/MBT_k4_minimum_enrichment.cov.RData"))
  h3ac.3h <- get_load(paste0(rdata_base_path(), "enrichment/MBT_h3ac_minimum_enrichment.cov.RData"))

  genes.list  <- list(sw.dev.low=genes.sw.dev.lo, sw.dev.high=genes.sw.dev.hi, sw.mat=genes.sw.mat, fw=genes.fw)

  pol.samples  <- list(tp2h=pol.2h, tp3h=pol.3h)
  tbp.samples  <- list(tp2h=tbp.2h, tp3h=tbp.3h)
  k4.samples   <- list(tp2h=k4.2h,  tp3h=k4.3h)
  h3ac.samples <- list(tp3h=h3ac.3h)

  reads.pol  <- get_enrichment_matrix_list(pol.samples,  genes.list)
  reads.tbp  <- get_enrichment_matrix_list(tbp.samples,  genes.list)
  reads.k4   <- get_enrichment_matrix_list(k4.samples,   genes.list)
  reads.h3ac <- get_enrichment_matrix_list(h3ac.samples, genes.list)
  
  varnames <- c("reads.pol", "reads.tbp", "reads.k4", "reads.h3ac")
  save(list=varnames, file=cache_file)
}

# calculate maximum enrichment (99th percentile)
pol.normalization  <- find_upper_threshold_for_samples(reads.pol,  quantile.threshold=0.99)
tbp.normalization  <- find_upper_threshold_for_samples(reads.tbp,  quantile.threshold=0.99)
k4.normalization   <- find_upper_threshold_for_samples(reads.k4,   quantile.threshold=0.99)
h3ac.normalization <- find_upper_threshold_for_samples(reads.h3ac, quantile.threshold=0.99)

# K4 normalization at 2h should use the quantile value for 3h (there is only noise at 2h)
k4.normalization$tp2h <- k4.normalization$tp3h

# scale so that 1 == maximum enrichment, 0 == background
reads.pol.n  <- normalize_matrices_by_threshold(reads.pol,  pol.normalization)
reads.tbp.n  <- normalize_matrices_by_threshold(reads.tbp,  tbp.normalization)
reads.k4.n   <- normalize_matrices_by_threshold(reads.k4,   k4.normalization)
reads.h3ac.n <- normalize_matrices_by_threshold(reads.h3ac, h3ac.normalization)

# randomly order the 4 gene lists
# for reproducibility, set random seed
set.seed(12345)

fw.order          <- sample(1:nrow(reads.pol.n$tp2$fw),          size=nrow(reads.pol.n$tp2$fw))
sw.mat.order      <- sample(1:nrow(reads.pol.n$tp2$sw.mat),      size=nrow(reads.pol.n$tp2$sw.mat))
sw.dev.high.order <- sample(1:nrow(reads.pol.n$tp2$sw.dev.high), size=nrow(reads.pol.n$tp2$sw.dev.high))
sw.dev.low.order  <- sample(1:nrow(reads.pol.n$tp2$sw.dev.low),  size=nrow(reads.pol.n$tp2$sw.dev.low))

gene.order <- list(fw=fw.order, sw.dev.high=sw.dev.high.order, sw.dev.low=sw.dev.low.order, sw.mat=sw.mat.order)

reads.pol.n  <- reorder_genelist(reads.pol.n,  gene.order)
reads.tbp.n  <- reorder_genelist(reads.tbp.n,  gene.order)
reads.k4.n   <- reorder_genelist(reads.k4.n,   gene.order)
reads.h3ac.n <- reorder_genelist(reads.h3ac.n, gene.order)

# number of blanks to insert between each gene group
blank.count <- nrow(reads.pol.n$tp2$fw)

pol.all  <- combine_genelists(reads.pol.n,  blank.count)
tbp.all  <- combine_genelists(reads.tbp.n,  blank.count)
k4.all   <- combine_genelists(reads.k4.n,   blank.count)
h3ac.all <- combine_genelists(reads.h3ac.n, blank.count)

g1.legend <- generate_plot(pol.all$tp2, "pre-MBT Pol (4)", show.legend=TRUE)

g1 <- generate_plot(pol.all$tp2, "pre-MBT Pol (4)", show.legend=FALSE)
g2 <- generate_plot(tbp.all$tp2, "pre-MBT TBP (2)", show.legend=FALSE)
g3 <- generate_plot(k4.all$tp2,  "pre-MBT K4 (2)",  show.legend=FALSE)

g4 <- generate_plot(pol.all$tp3,  "MBT Pol (3)",  show.legend=FALSE)
g5 <- generate_plot(tbp.all$tp3,  "MBT TBP (2)",  show.legend=FALSE)
g6 <- generate_plot(k4.all$tp3,   "MBT K4 (2)",   show.legend=FALSE)
g7 <- generate_plot(h3ac.all$tp3, "MBT H3Ac (2)", show.legend=FALSE)

pdf("figure_1c_heatmap.pdf", width=4, height=9, onefile=T)

print(g1.legend)
print(g1)
print(g2)
print(g3)
print(g4)
print(g5)
print(g6)
print(g7)

dev.off()









