library(Biostrings)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(data.table)

source("shared_code/flybase.r")

genes <- flybase_with_custom_txs()
valid.chrs <- c("chr2L", "chr2LHet", "chr2R", "chr2RHet", "chr3L", "chr3LHet", "chr3R", "chr3RHet", "chr4", "chrX", "chrXHet", "chrYHet")
genes <- subset(genes, chr %in% valid.chrs)

piv <- which(genes$strand == 1)
niv <- which(genes$strand == -1)

results.df <- NULL

zelda.motif <- "YAGGTAR"

message("Motif: ", zelda.motif) 
  
motif     <- DNAString(zelda.motif)
new_start <- -2000
new_end   <- 0

# determine promoter region to scan
genes$pstart <- NA
genes$pend   <- NA

genes[piv, ] <- transform(genes[piv, ], pstart = start + new_start,
                                        pend   = start + new_end)

genes[niv, ] <- transform(genes[niv, ], pstart = end - new_end,
                                        pend   = end - new_start)

chr_limits <- seqlengths(Dmelanogaster)

igenes <- subset(genes, pstart > 0 & pstart <= chr_limits[as.character(chr)] & 
                        pend   > 0 & pend   <= chr_limits[as.character(chr)])

message(" - fetching sequences")
forward_seq <- getSeq(Dmelanogaster, names=igenes$chr, start=igenes$pstart, end=igenes$pend, strand="+", as.character=FALSE)
reverse_seq <- getSeq(Dmelanogaster, names=igenes$chr, start=igenes$pstart, end=igenes$pend, strand="-", as.character=FALSE)

names(forward_seq) <- paste(igenes$fb_tx_id, "_f", sep="")
names(reverse_seq) <- paste(igenes$fb_tx_id, "_r", sep="")

promoter_seqs <- c(forward_seq, reverse_seq)

message(" - finding motif matches")
mindex_0mm <- vmatchPattern(motif, promoter_seqs, fixed=FALSE, max.mismatch=0)
mindex_1mm <- vmatchPattern(motif, promoter_seqs, fixed=FALSE, max.mismatch=1)

df.0mm <- data.frame(mm=0, fb_tx_id=names(mindex_0mm), count=countIndex(mindex_0mm))
df.1mm <- data.frame(mm=1, fb_tx_id=names(mindex_1mm), count=countIndex(mindex_1mm))

message("Combining matches on forward and reverse strand...")
motif.df <- rbind(df.0mm, df.1mm)
motif.df$fb_tx_id <- gsub("_.$", "", motif.df$fb_tx_id)

motif.dt <- data.table(motif.df)
motif.dt <- motif.dt[ , list(count=sum(count)), by=list(mm, fb_tx_id)]
motif.df <- as.data.frame(motif.dt)
motif.df$motif_name <- "zelda"

zelda.counts <- motif.df
save(zelda.counts, file="zelda_motif_counts.RData")

