library(data.table)

message("Reading in exon/intron mRNA lines...")
gff <- read.delim(pipe("gzcat dmel-all-no-analysis-r5.47.gff.gz | awk '{ if ($3 == \"exon\" || $3 == \"intron\") print }' | grep 'parent_type=mRNA'"), stringsAsFactors=F, header=F)

gff$fb_tx_id <- gsub("^.*;Parent=(.*?);.*$", "\\1", gff$V9)

gff.new <- transform(gff[rep(seq(nrow(gff)), sapply(v.s <- strsplit(gff$fb_tx_id, split=","), length)),], fb_tx_id=unlist(v.s))
gff.new$V9 <- NULL

gff.new <- transform(gff.new, width = V5 - V4 + 1)
names(gff.new)[which(names(gff.new) == "V3")] <- "feature"

gff.dt <- data.table(gff.new[, c("feature", "width", "fb_tx_id")])

counts <- as.data.frame(gff.dt[, list(exon_sum=sum(width[feature == "exon"]), intron_sum=sum(width[feature == "intron"])), by=fb_tx_id])

txs <- get(load("fb.transcripts.r5.47.RData"))

txs <- merge(txs, counts)
txs <- transform(txs, total_width = end - start + 1)

intron_exon_counts.df <- txs
save(intron_exon_counts.df, file="intron_exon_counts.df.RData")
