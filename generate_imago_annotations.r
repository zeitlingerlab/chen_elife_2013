
insitu <- read.csv("imago/insitu_annot.csv", stringsAsFactors=F, header=F)

names(insitu) <- c("fb_symbol", "fb_cg_id", "fb_gene_id", "stage_id", "go_term")

stages <- data.frame(stringsAsFactors=F, stage_id=1:6, stage_name=c("stage01-03", "stage04-06", "stage07-08", "stage09-10", "stage11-12", "stage13-16"))
insitu <- merge(insitu, stages)

source("shared_code/flybase.r")

flybase.genes <- unique(flybase_with_custom_txs()[, c("fb_gene_id", "fb_symbol", "fb_cg_id")])

insitu$fb_gene_id[insitu$fb_gene_id == ""] <- NA

# match symbol to fb_gene_id if fb_gene_id is missing

missing_id <- which(is.na(insitu$fb_gene_id))
message("Missing fb_gene_id: ", length(missing_id))
message("Matching using fb_symbol...")
insitu$fb_gene_id[missing_id] <- flybase.genes$fb_gene_id[match(insitu$fb_symbol[missing_id], flybase.genes$fb_symbol)]

missing_id <- which(is.na(insitu$fb_gene_id))
message("Still missing fb_gene_id: ", length(missing_id))
message("Matching using fb_cg_id...")
insitu$fb_gene_id[missing_id] <- flybase.genes$fb_gene_id[match(insitu$fb_cg_id[missing_id], flybase.genes$fb_cg_id)]

missing_id <- which(is.na(insitu$fb_gene_id))
message("Still missing fb_gene_id: ", length(missing_id))

message("Matching with fb_cg_id as fb_symbol...")
insitu$fb_gene_id[missing_id] <- flybase.genes$fb_gene_id[match(insitu$fb_cg_id[missing_id], flybase.genes$fb_symbol)]

missing_id <- which(is.na(insitu$fb_gene_id))
message("Still missing fb_gene_id: ", length(missing_id))

insitu$stage_id <- NULL
save(insitu, file="imago/insitu.RData")
