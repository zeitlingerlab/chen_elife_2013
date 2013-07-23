library(GenomicRanges)
library(lattice)

# given a data frame of intervals (columns chr, region_start and region_end), returns a matrix
# of values in those intervals with row names taken from id.column
read_matrix <- function(loc.df, cov.object, id.column="fb_gene_id") {
  if(length(unique(loc.df$region_end - loc.df$region_start)) > 1) stop("Provided intervals are not all the same length")
  heatmap.rd         <- RangedData(IRanges(start=loc.df$region_start, end=loc.df$region_end), space=loc.df$chr, id=loc.df[, id.column])
  heatmap.view       <- RleViewsList(rleList=cov.object[names(heatmap.rd)], rangesList=ranges(heatmap.rd))
  heatmap.reads      <- viewApply(heatmap.view, as.numeric)
  heatmap.matrix     <- matrix(unlist(sapply(heatmap.reads, sapply, as.numeric)), nrow=nrow(loc.df), byrow=TRUE)
  rownames(heatmap.matrix) <- as.character(as.data.frame(heatmap.rd)$id)
  heatmap.matrix
}

# given a data frame listing genes (columns chr, start, end and strand), returns a matrix
# of read values for an equal-sized region in each gene defined by upstream and downstream
get_enrichment_matrix <- function(cov.object, genes, id.column="fb_gene_id", upstream=200, downstream=800) {
  if(!identical(sort(unique(genes$strand)), c(-1, 1))) stop("genes data frame should have a strand column with 1 and -1")
  chr_lengths <- sapply(cov.object, length)

  original.genes.count <- nrow(genes)
  genes <- transform(genes, region_start = ifelse(strand == 1, start - upstream, end - downstream + 1),
                            region_end   = ifelse(strand == 1, start + downstream - 1, end + upstream))

  genes <- subset(genes, region_start > 0 & region_end <= chr_lengths[as.character(chr)])
  
  if(nrow(genes) != original.genes.count) warning(original.genes.count - nrow(genes), " gene(s) were removed due to chromosome boundary")
  
  niv <- which(genes$strand == -1)
  piv <- which(genes$strand == 1)
  
  reads.p <- NULL
  reads.n <- NULL
  
  if(length(piv) > 0) {
    reads.p <- read_matrix(genes[piv, ], cov.object, id.column)
  }

  if(length(niv) > 0) {
    reads.n <- read_matrix(genes[niv, ], cov.object, id.column)
    reads.n <- reads.n[, ncol(reads.n):1]
  }

  reads <- rbind(reads.p, reads.n)
  reads
}

# given a list of samples (coverage objects) and a list of gene data frames, generate
# read matrices for every combination 
get_enrichment_matrix_list <- function(samples, genes, ...) {
  sample.names <- names(samples)
  sample.counter <- 1
  genelist.names <- names(genes)
  lapply(samples, function(sample) {
    # for each coverage object (sample), loop through gene lists
    message("get_enrichment_matrix_list() On sample: ", sample.names[sample.counter])
    sample.counter <<- sample.counter + 1
    genelist.counter <- 1
    lapply(genes, function(genelist) {
      message(" - gene list: ", genelist.names[genelist.counter])
      genelist.counter <<- genelist.counter + 1
      get_enrichment_matrix(sample, genelist, ...)
    })
  })
}


normalize_matrix <- function(reads, value.limit) {
  reads[reads < 0] <- 0
  reads <- reads / value.limit
  reads[reads > 1] <- 1
  reads
}

generate_plot <- function(reads, plot.title, show.legend=TRUE) {
  rgb.palette <- colorRampPalette(c("white", "blue", "red"), space = "rgb")
  message("plot: ", plot.title)
  rownames(reads) <- NULL
  levelplot(t(reads), main=plot.title, xlab="", ylab="", col.regions=rgb.palette(120), useRaster=TRUE, cuts=16, colorkey=show.legend)
}

# For each sample in sample.list, combine all the matrices for that sample and calculate a quantile
find_upper_threshold_for_samples <- function(sample.list, quantile.threshold=0.99) {
  lapply(sample.list, function(s) quantile(c(do.call(rbind, s)), quantile.threshold, na.rm=T))
}

# apply normalization to all samples in sample.list using the matching threshold in threshold.list
normalize_matrices_by_threshold <- function(sample.list, threshold.list) {
  for(s in names(sample.list)) {
    message("apply_normalization_targets() sample: ", s)
    for(n in names(sample.list[[s]])) {
      value.target <- threshold.list[[s]]
      message(" - target for ", n, ": ", value.target)
      sample.list[[s]][[n]] <- normalize_matrix(sample.list[[s]][[n]], value.target)
    }
  }
  sample.list
}

# re-orders each matrix in sample.list using the matching order in order.list
reorder_genelist <- function(sample.list, order.list) {
  for(s in names(sample.list)) {
    for(g in names(sample.list[[s]])) {
      if(g %in% names(order.list)) {
        message("Re-ordering ", g, " in ", s)
        sample.list[[s]][[g]] <- sample.list[[s]][[g]][order.list[[g]], ]
        if(nrow(sample.list[[s]][[g]]) != length(order.list[[g]])) stop("reorder_genelist(): provided order does not have enough values for specified gene list")
      }
    }
  }
  sample.list
}

# combine a list of matrices, but add blanks between them
rbind_matrices_with_blanks <- function(mlist, blanks) {
  first_one <- TRUE
  m <- NULL
  for(i in 1:length(mlist)) {
    if(first_one == TRUE) {
      m <- mlist[[i]]
      first_one <- FALSE
    } else {
      m <- rbind(m, blanks, mlist[[i]])
    }
  }
  m
}

# merge the gene matrices for each sample
combine_genelists <- function(sample.list, empty.rows.between) {
  m.columns <- ncol(sample.list[[1]][[1]])
  blank.m <- matrix(rep(NA, times=m.columns*empty.rows.between), ncol=m.columns)
  
  lapply(sample.list, function(mlist) rbind_matrices_with_blanks(mlist, blank.m))    
}

# reverse the column order of a list of matrices
reverse_list_of_matrices <- function(mlist) {
  lapply(mlist, function(m) m[, ncol(m):1])
}
