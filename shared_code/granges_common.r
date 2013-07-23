library(GenomicRanges)

regionSums <- function(regions, cvg) {
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewSums(
             Views(cvg, as(regions, "RangesList"))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionWhichMaxs <- function(regions, cvg) {
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewWhichMaxs(
             Views(cvg, as(regions, "RangesList"))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

total_signal <- function(cov) {
  sum(as.numeric(sapply(cov, function(x) sum(as.numeric(x)))))
}

get_load <- function(filename) {
  message("Loading ", filename, " ... ", appendLF=FALSE)
  o <- updateObject(get(load(filename)))
  message("OK")
  o
}

