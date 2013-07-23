
ids <- read.delim("fbgn_annotation_ID_fb_2012_05.tsv.gz", stringsAsFactors=F, header=F, skip=5)
names(ids) <- c("fb_symbol", "fb_gene_id", "prev_id", "cg_id", "other")
ids$fb_symbol <- NULL
ids$cg_id <- NULL
ids$other <- NULL

ids.expand <- subset(ids, prev_id != "")
ids.expand <- transform(ids.expand[rep(seq(nrow(ids.expand)), sapply(v.s <- strsplit(ids.expand$prev_id, split=","), length)),], prev_id=unlist(v.s))

ids.same <- data.frame(stringsAsFactors=F, fb_gene_id = unique(as.character(ids$fb_gene_id)), 
                                           prev_id    = unique(as.character(ids$fb_gene_id)))

ids.all <- rbind(ids.same, ids.expand)

fbidmap <- ids.all
save(fbidmap, file="fbidmap.RData")
