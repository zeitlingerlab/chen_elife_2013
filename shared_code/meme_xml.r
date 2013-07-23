library(XML)
library(seqLogo)

extract_pwm_from_motif_node <- function(motif_node) {
  alpha.m <- motif_node[["probabilities"]][["alphabet_matrix"]]
  alpha.m <- as.numeric(xmlSApply(alpha.m, xmlSApply, xmlValue))
  alpha.m <- matrix(alpha.m, nrow=4, byrow=FALSE)
  # assume ACGT order (MEME XML)
  rownames(alpha.m) <- c("A", "C", "G", "T")
  alpha.m
}

extract_motif_attrs_from_motif_node <- function(motif_node) {
  attrs <- xmlAttrs(motif_node)
  as.list(attrs)
}

extract_contributing_sites_table <- function(motif_node) {
  sites.m <- xmlSApply(motif_node[["contributing_sites"]], xmlAttrs)
  row.seq    <- which(rownames(sites.m) == "sequence_id")
  row.pos    <- which(rownames(sites.m) == "position")
  row.strand <- which(rownames(sites.m) == "strand")
  data.frame(stringsAsFactors=F, sequence_id = sites.m[row.seq, ],
                                 position = as.integer(sites.m[row.pos, ]),
                                 strand   = ifelse(as.character(sites.m[row.strand, ]) == "plus", "+", "-"))
}

extract_sequence_name_map <- function(training_set_node) {
  ts.children <- xmlChildren(training_set_node)
  sequence.ids <- which(names(ts.children) == "sequence")
  seq.m <- do.call(rbind, lapply(ts.children[sequence.ids], xmlAttrs))
  rownames(seq.m) <- NULL
  seq.m <- as.data.frame(stringsAsFactors=F, seq.m)[, c("id", "name")]
  names(seq.m) <- c("sequence_id", "name")
  seq.m
}

parse_meme_xml <- function(filename) {
  doc <- xmlRoot(xmlTreeParse(filename))
  motifs <- xmlSApply(doc[["motifs"]], extract_pwm_from_motif_node)
  names(motifs) <- paste("motif", as.character(1:length(motifs)), sep="_")  

  motif_info <- xmlSApply(doc[["motifs"]], extract_motif_attrs_from_motif_node)
  name.row <- which(rownames(motif_info) == "id")
  colnames(motif_info) <- motif_info[name.row, ]

  motif_sites <- xmlApply(doc[["motifs"]], extract_contributing_sites_table)
  names(motif_sites) <- paste("motif", as.character(1:length(motifs)), sep="_")  
    
  input_sequence_count <- length(which(xmlSApply(doc[["training_set"]], xmlName) == "sequence"))
  
  sequence_name_map <- extract_sequence_name_map(doc[["training_set"]])
  
  list(motifs            = motifs, 
       motif_info        = as.data.frame(motif_info),
       motif_sites       = motif_sites,
       sequence_count    = input_sequence_count,
       sequence_name_map = sequence_name_map)
}

motif_sites <- function(meme, motif_name) {
  name.map <- meme$sequence_name_map
  sites    <- meme$motif_sites[[motif_name]]
  merge(name.map, sites)
}

