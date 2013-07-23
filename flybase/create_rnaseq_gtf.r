
message("Reading in exon mRNA lines...")
exons <- read.delim(pipe("gzcat dmel-all-no-analysis-r5.47.gff.gz | awk '{ if ($3 == \"exon\") print }' | grep 'parent_type=mRNA'"), stringsAsFactors=F, header=F)

exons$gene_ids       <- gsub("^.*ID=(.*):.*;Name=.*$",         "\\1", exons$V9)
exons$transcript_ids <- gsub("^.*Parent=(.*);parent_type=.*$", "\\1", exons$V9)

# some gene_ids are CG ids instead of FB ids
message("Reading annotation ID map file...")
fbnames <- read.delim("fbgn_annotation_ID_fb_2012_05.tsv.gz", stringsAsFactors=F, header=F, skip=5)
names(fbnames) <- c("symbol", "fb_id", "secondary_fb_id", "cg_id", "other")

cgids <- grep("^CG", exons$gene_ids)
exons$gene_ids[cgids] <- fbnames$fb_id[match(exons$gene_ids[cgids], fbnames$cg_id)]

message("Expanding exon lines...")
exons <- transform(exons[rep(seq(nrow(exons)), sapply(v.s <- strsplit(exons$transcript_ids, split=","), length)),], flybase_transcript_id=unlist(v.s))
exons <- transform(exons, V9 = paste("gene_id \"", gene_ids, "\"; transcript_id \"", flybase_transcript_id, "\"", sep=""))

#spikes <- read.delim("spikes.gtf", stringsAsFactors=F, header=F)

message("Writing file...")
write.table(exons[, 1:9], file="fb-r5.47.gtf", row.names=F, col.names=F, sep="\t", quote=F)

# check with cufflink's gffread to remove bad transcripts
message("Checking with gffread...")
dups <- system("~/apps/cufflinks2/gffread fb-r5.47.gtf 2>&1", intern=T)

tx_dups <- unique(gsub("^GFF Error at (.*) ...: .*$", "\\1", dups))
message("Removing ", length(tx_dups), " transcript_id(s)")
message(paste(tx_dups, collapse=", "))

exons <- subset(exons, !(flybase_transcript_id %in% tx_dups))
message("Re-writing file...")
write.table(exons[, 1:9], file="fb-r5.47.clean.gtf", row.names=F, col.names=F, sep="\t", quote=F)
