suppressPackageStartupMessages(library(GenomicRanges))

args <- commandArgs(trailingOnly=TRUE)

samples <- list()

get_load <- function(f) {
  message("Loading: ", f)
  updateObject(get(load(f)))
}

for(f in args[-length(args)]) {
  s.cov <- get_load(f)
  samples <- c(samples, list(s.cov))
}

min.cov <- do.call(pmin, samples)

save(min.cov, file=args[length(args)])
