options(scipen=999)

#genes <- get(load("fb.transcripts.r5.46.custom.RData"))
#
#genes$strand <- ifelse(genes$strand == 1, "+", "-")
#genes$score <- 0
## chr start stop name score strand
#
#write.table(genes[, c("chr", "start", "end", "fb_symbol", "score", "strand")], file="transcripts_igv_custom.bed", quote=F, row.names=F, col.names=F, sep="\t")


genes <- get(load("fb.transcripts.r5.47.RData"))

genes$strand <- ifelse(genes$strand == 1, "+", "-")
genes$score <- 0
# chr start stop name score strand

write.table(genes[, c("chr", "start", "end", "fb_symbol", "score", "strand")], file="transcripts_igv_flybase.bed", quote=F, row.names=F, col.names=F, sep="\t")

