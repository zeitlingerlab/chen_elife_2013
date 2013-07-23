suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-t", "--treatment"), 
              type="character",
              default=NA,
              help="Treatment coverage file"),
  make_option(c("-c", "--control"), 
              type="character",
              default=NA,
              help="Control coverage file"),
  make_option(c("-n", "--name"),
              type="character",
              default=NA,
              help="Name of resulting coverage object"),
  make_option(c("-b", "--bigwig"),
              type="character",
              default=NA,
              help="Name of resulting BigWig file"),
  make_option(c("-w", "--window"),
              type="integer",
              default=201,
              help="Sliding window size"),
  make_option(c("-s", "--skipbigwig"),
              action="store_true",
              default=FALSE,
              help="Skip creation of BigWig file"),
  make_option(c("-l", "--transform"),
              type="character",
              default="log2",
              help="Transform to apply to enrichment values (log2 or linear)")
  )

opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$treatment)) {
  message("No treatment file specified.")
  q(status=1)
}

if(is.na(opt$control)) {
  message("No control file specified.")
  q(status=1)
}

if(is.na(opt$name)) {
  message("No output name specified.")
  q(status=1)
}

if(! opt$transform %in% c("log2", "linear")) {
  message("Transform type must be either 'log2' or 'linear'")
  q(status=1)
}

suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F))

load_cov <- function(filename) {
  message("Loading: ", filename)
  get(load(filename))
}

ip.cov <- updateObject(load_cov(opt$treatment))
bg.cov <- updateObject(load_cov(opt$control))

common.chrs <- intersect(names(ip.cov), names(bg.cov))
ip.cov <- ip.cov[names(ip.cov) %in% common.chrs]
bg.cov <- bg.cov[names(bg.cov) %in% common.chrs]

bw.target <- ifelse(is.na(opt$bigwig), paste(opt$name, ".bw", sep=""), opt$bigwig)

convert_to_running_window_sums <- function(cov, window_size, add_background=0) {
  for(chr in names(cov)) {
    #message("chr: ", chr)
    if(add_background != 0) cov[[chr]] <- cov[[chr]] + add_background
    cov[[chr]] <- runsum(Rle(as.numeric(cov[[chr]])), window_size, endrule="constant")
  }
  cov
}

convert_to_enrichment <- function(cov.ip, cov.bg, window_size=101, add_background=1, min_enrichment=0, transform_function=identity) {
  ip.sum <- sum(as.numeric(sapply(cov.ip, function(x) sum(as.numeric(x)))))
  bg.sum <- sum(as.numeric(sapply(cov.bg, function(x) sum(as.numeric(x)))))
  
  message("Processing treatment...")
  ip.rw <- convert_to_running_window_sums(cov.ip, window_size)
  message("Processing control...")
  bg.rw <- convert_to_running_window_sums(cov.bg, window_size, add_background)
  
  message("Calculating enrichments...")
  for(chr in names(ip.rw)) {
    ip.rw[[chr]] <- Rle(pmax(round(as.numeric(transform_function((as.numeric(ip.rw[[chr]]) / as.numeric(ip.sum)) / (as.numeric(bg.rw[[chr]]) / as.numeric(bg.sum)))), 2), min_enrichment))
  }
  ip.rw
}

cov_to_bigwig <- function(cov, filename) {
  cov.lengths <- sapply(cov, length)
  cov.ranged  <- as(cov, "RangedData")
  seqlengths(cov.ranged) <- cov.lengths
 
  export(cov.ranged, filename, dataFormat="bedGraph")
}

transform_function <- get(ifelse(opt$transform == "log2", "log2", "identity"))

e.cov <- convert_to_enrichment(cov.ip=ip.cov, cov.bg=bg.cov, window_size=opt$window, transform_function=transform_function)
rm(ip.cov)
rm(bg.cov)
nothing <- gc()

assign(opt$name, e.cov)

message("Saving enrichments...")
save(list=opt$name, file=paste(opt$name, ".cov.RData", sep=""))

if(!opt$skipbigwig) {
  message("Saving bigwig...")
  nothing <- cov_to_bigwig(e.cov, bw.target)
}

