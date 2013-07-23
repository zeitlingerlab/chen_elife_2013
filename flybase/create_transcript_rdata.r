library(plyr)

message("Reading in RNA lines...")
tx <- read.delim(pipe("gzcat dmel-all-no-analysis-r5.47.gff.gz | awk '{ if ($3 == \"mRNA\" || $3 == \"ncRNA\" || $3 == \"rRNA\" || $3 == \"snRNA\" || $3 == \"snoRNA\") print }'"), stringsAsFactors=F, header=F)

tx$transcript_id <- gsub("^ID=(.*?);.*$",        "\\1", tx$V9)
tx$gene_id       <- gsub("^.*;Parent=(.*?);.*$", "\\1", tx$V9)

fb.transcripts <- data.frame(stringsAsFactors=F, fb_gene_id = tx$gene_id, 
                                                 fb_tx_id   = tx$transcript_id, 
                                                 chr        = paste("chr", tx$V1, sep=""),
                                                 start      = tx$V4,
                                                 end        = tx$V5,
                                                 strand     = ifelse(tx$V7 == "+", 1, -1),
                                                 type       = tx$V3)

# Remove "chrdmel_mitochondrion_genome"
fb.transcripts <- subset(fb.transcripts, chr != "chrdmel_mitochondrion_genome")

# select one transcript among those with the same coordinates
selected.tx <- ddply(fb.transcripts, .var=.(fb_gene_id, chr, start, end, strand, type), summarize, selected_transcript_id=fb_tx_id[1], .progress="text")

fb.transcripts <- subset(fb.transcripts, fb_tx_id %in% selected.tx$selected_transcript_id)

message("Reading annotation ID map file...")
fbnames <- read.delim("fbgn_annotation_ID_fb_2012_05.tsv.gz", stringsAsFactors=F, header=F, skip=5)[, c(1,2,4)]
names(fbnames) <- c("fb_symbol", "fb_gene_id", "fb_cg_id")

# pre_miRNA
pm <- read.delim(pipe("gzcat dmel-all-no-analysis-r5.47.gff.gz | awk '{ if ($3 == \"pre_miRNA\") print }'"), stringsAsFactors=F, header=F)
pm$transcript_id <- gsub("^ID=(.*?);.*$",        "\\1", pm$V9)
pm$gene_id       <- gsub("^.*;Parent=(.*?);.*$", "\\1", pm$V9)

pm.df <- data.frame(stringsAsFactors=F, fb_gene_id = pm$gene_id, 
                                        fb_tx_id   = pm$transcript_id, 
                                        chr        = paste("chr", pm$V1, sep=""),
                                        start      = pm$V4,
                                        end        = pm$V5,
                                        strand     = ifelse(pm$V7 == "+", 1, -1),
                                        type       = pm$V3)

pm.add <- data.frame(stringsAsFactors=F, fb_gene_id = "FBgn0262380",
                                         fb_tx_id   = "FBtr0304223",
                                         chr        = "chr2R",
                                         start      = 15548221,
                                         end        = 15549279,
                                         strand     = -1,
                                         type       = "pre_miRNA_cluster")
pm.df <- subset(pm.df, ! (chr == "chr2R" & start >= 15548221 & end <= 15549279))

fb.transcripts <- rbind(pm.df, pm.add, fb.transcripts)

fb.transcripts.r5.47 <- merge(fbnames, fb.transcripts, all.y=T)

save(fb.transcripts.r5.47, file="fb.transcripts.r5.47.RData")
