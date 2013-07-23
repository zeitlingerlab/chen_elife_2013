library(ggplot2)
library(scales)
library(gtools)
library(matrixStats)
library(fastcluster)
library(simpleaffy)

source("shared_code/stat_tests.r")

pe <- get(load("promoter_elements/pe.0mm.RData"))
all.fb.ids <- unique(as.character(pe$fb_tx_id))

add_list <- function(genes, name) {
  genes <- unique(as.character(genes))
  gene.lists <<- c(gene.lists, list(genes))
  #names(gene.lists)[length(gene.lists)] <- paste(name, "\n", prettyNum(length(genes), big.mark=","), " genes", sep="")
  names(gene.lists)[length(gene.lists)] <- name
  gene.lists <<- gene.lists
  return(TRUE)
}

pe_enrichments <- function(pe, genes, name) {

  results <- NULL

  for(i in 2:ncol(pe)) {
    ftest <- fisher_test_2x2(genes, pe$fb_tx_id[pe[, i] == TRUE], pe$fb_tx_id)
    ftest$element <- names(pe)[i]
    ftest$gene_group <- name
    ftest$gene_group_size <- length(genes)
    results <- rbind(results, ftest)
  }
  results
}

pe_heatmap_for_tx_lists <- function(lists, name, cluster.groups=TRUE, manual.group.order=c()) {
  
  gene.lists <<- NULL
  
  message("Creating gene lists...")
  for(listname in names(lists)[-1]) {
    genes <- lists[which(lists[, listname] == TRUE), ]$fb_tx_id
    message(listname, ": ", length(genes), " gene(s)")
    add_list(genes, listname)
  }

  plot.data <- NULL
  for(i in 1:length(gene.lists)) {
    list.name <- names(gene.lists)[i]
    list.genes <- gene.lists[[i]]
    message(list.name, " (", length(list.genes), " genes)")
    pe.e <- pe_enrichments(pe, list.genes, list.name)
    plot.data <- rbind(plot.data, pe.e)
  }
  plot.data$gene_group <- factor(plot.data$gene_group, levels=names(gene.lists))

  plot.data <- transform(plot.data, significant = ifelse(pvalue < 0.05, "*", ""))
  plot.data <- transform(plot.data, pvalue = -1 * log10(pvalue))
  plot.data <- transform(plot.data, pvalue = ifelse(test_type == "Enrichment", pvalue, -1 * pvalue))

  # cap pvalues at 25
  plot.data.original <- plot.data
  plot.data <- transform(plot.data, pvalue = ifelse(pvalue > 25, 25, pvalue))
  plot.data <- transform(plot.data, pvalue = ifelse(pvalue < -25, -25, pvalue))

  # cluster:

  # for clustering

  data.wide <- reshape(data=plot.data.original[, c("gene_group", "pvalue", "element")], v.names="pvalue", idvar="gene_group", timevar="element", direction="wide")
  row.names(data.wide) <- data.wide$gene_group
  data.wide <- data.wide[, -1]
  #write.table(data.wide, file="pe_pvalues_cluster.txt", quote=F, sep="\t", row.names=T, col.names=NA)

  #group.order   <- row.names(data.wide)[order.dendrogram(as.dendrogram(hclust(dist(as.matrix(data.wide)))))]
  #element.order <- gsub("pvalue\\.", "", row.names(t(data.wide))[order.dendrogram(as.dendrogram(hclust(dist(t(as.matrix(data.wide))))))])

  colnames(data.wide) <- gsub("pvalue\\.", "", colnames(data.wide))

  # pearson
  col.d <- as.dendrogram(standard.pearson(as.matrix(data.wide)))
  row.d <- as.dendrogram(standard.pearson(t(as.matrix(data.wide))))

  #col.d <- as.dendrogram(hclust(dist(t(as.matrix(data.wide)), method="manhattan"), method="complete"))
  #row.d <- as.dendrogram(hclust(dist(as.matrix(data.wide), method="manhattan"), method="complete"))

  element.order <- row.names(t(as.matrix(data.wide)))[order.dendrogram(col.d)]
  group.order   <- row.names(data.wide)[order.dendrogram(row.d)]


  if(cluster.groups == TRUE) plot.data$gene_group <- factor(plot.data$gene_group, levels=rev(group.order))
  if(length(manual.group.order) > 0) plot.data$gene_group <- factor(plot.data$gene_group, levels=rev(manual.group.order))

  elements <- c("TATA", "Inr", "DPE", "PB", "MTE", "DRE", "Motif1", "Motif6", "Motif7", "GAGA", "Zelda")
  plot.data$element <- factor(plot.data$element, levels=elements)
  plot.data <- subset(plot.data, !is.na(element))
    
  plot.data <- transform(plot.data, enrichment = ifelse(test_type == "Depletion", -1 * (1 / enrichment), enrichment))
  plot.data$enrichment <- pmax(-4, plot.data$enrichment)
  plot.data$enrichment <- pmin(4, plot.data$enrichment)

  e.limit.up   <- max(plot.data$enrichment)
  e.limit.down <- min(plot.data$enrichment)

  plot.data <- transform(plot.data, gene_group_label = paste0(gene_group, "\n", prettyNum(gene_group_size, big.mark=",", preserve.width="individual"), " genes"))
  plot.data <- plot.data[order(plot.data$gene_group), ]
  new_labels <- unique(plot.data$gene_group_label)
  plot.data$gene_group <- factor(plot.data$gene_group, labels=new_labels)

  e.limit.up <- max(plot.data$enrichment)
  e.limit.down <- min(plot.data$enrichment)

  g.e <- ggplot(plot.data, aes(x=element, y=gene_group, fill=enrichment)) + 
          geom_tile() + 
          geom_text(aes(label=significant), color="white") +
          theme_bw() +
          coord_equal() +
          #scale_fill_gradient2("Enrichment", space="rgb", low="black", mid="white", high="#FC8F00", limits=c(-1*e.limit, e.limit), midpoint=0, guide=guide_colorbar()) +
          scale_fill_gradientn(name="Enrichment", space="rgb", 
                               values=c(e.limit.down, -0.5, 0.5, e.limit.up), 
                               colours=c("#000000", "#cccccc", "#cccccc", "#FC8F00"), 
                               rescaler=function(x,...) x, oob=identity,
                               limits=c(e.limit.down, e.limit.up), guide=guide_colorbar()) +
          labs(x="", y="", title=name) +
          theme(panel.grid.minor=element_blank(),
                panel.grid.major=element_blank()) 

  list("plot"=g.e, "data"=plot.data)
}



