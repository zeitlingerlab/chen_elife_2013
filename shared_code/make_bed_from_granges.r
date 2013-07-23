library(GenomicRanges)
library(rtracklayer)
library(xlsx)
library(stringr)

get_load <- function(filename) {
  message("Loading ", filename, " ... ", appendLF=FALSE)
  o <- updateObject(get(load(filename)))
  message("OK")
  o
}

samples.df <- read.xlsx("/Volumes/ZeitlingerLab/Projects/early\ embryos/Pol\ II/GEO\ Submission/samples.xlsx", sheetIndex=1, stringsAsFactors=F)
samples.df$Rdata <- str_trim(samples.df$Rdata)

for(i in 1:nrow(samples.df)) {
  ranges.file <- paste("~/rdata/chromatin_early_embryo/", samples.df$Rdata[i], ".ranges.RData", sep="")
  gr <- get_load(ranges.file)
  bed.name <- gsub("\\.fastq\\.gz$", ".bed", samples.df$Sequence.file[i])
  message("Writing: ", bed.name)
  export(gr, paste("bed/", bed.name, sep=""))
}


