library(GenomicRanges)
library(matrixStats)

get_load <- function(filename) {
  message("Loading ", filename, " ... ", appendLF=FALSE)
  o <- updateObject(get(load(filename)))
  message("OK")
  o
}

calc_total_signal_for_coverage_object <- function(cov.object) {
  sum(as.numeric(sapply(cov.object, function(x) sum(as.numeric(x)))))
}

check_gene_list_for_proper_columns <- function(gene_list) {
  if(sum(c("chr", "start", "stop", "strand") %in% names(gene_list)) != 4) {
    stop("gene_list is missing required columns (chr, start, stop, strand)")
  }
  uniq_strand_values <- unique(gene_list$strand)
  if(!identical(uniq_strand_values, uniq_strand_values[uniq_strand_values %in% c(-1, 1)])) {
    stop("gene_list$strand must be 1 or -1. Values found: ", paste(sort(unique(gene_list$strand)), sep=" "))
  }  
}

get_metagene_scaled <- function(sample.cov, gene_list, width, smooth=NULL, sample=NULL) {
  sample_object_name <- deparse(substitute(sample.cov))

  # update the provided object in case it was built with a previous version of IRanges
  sample.cov <- updateObject(sample.cov)
  
  message("Scaled metagene call for ", nrow(gene_list), " genes (", sample_object_name, ")")
  
  check_gene_list_for_proper_columns(gene_list)
    
  # extract reads for each gene and scale to width
  reads <- matrix(nrow=nrow(gene_list), ncol=width)
  for(i in 1:nrow(gene_list)) {
    if(gene_list$strand[i] == 1) {
      read_start <- gene_list$start[i]
      read_stop  <- gene_list$stop[i]
    } else if (gene_list$strand[i] == -1) {
      read_start <- gene_list$stop[i]
      read_stop  <- gene_list$start[i]
    }
    
    chr_max <- length(sample.cov[[as.character(gene_list$chr[i])]])
    if(read_start <= 0 | read_stop <= 0 | read_start > chr_max | read_stop > chr_max) {
      message("Coordinates exceed chromosome limits: index ", i, ", ", gene_list$chr[i], ", skipping.")
      next
    }
    this_gene <- as.numeric(sample.cov[[as.character(gene_list$chr[i])]][read_start:read_stop])
    if(length(which(is.na(this_gene))) > 0) {
      stop("NA values found: ", gene_list$chr[i], " ", read_start, "-", read_stop)
    }
    if(!is.null(smooth)) {
      #message("Smoothing: ", i, " (", length(this_gene), " values) ", appendLF=F)
      this_gene <- as.numeric(runmean(Rle(this_gene), smooth, "constant"))
      #message("done.")
    }
    this_gene.scaled <- approx(this_gene, n=width)$y
    reads[i, ] <- this_gene.scaled
  }
  results <- data.frame(position=1:width, avg_reads=colMeans(reads))
  if(!is.null(sample)) { 
    results$sample <- sample
  }
  results
}

pad_with_nas <- function(x, target_length, reverse=FALSE) {
  if(reverse == TRUE) x <- rev(x)
  if(length(x) < target_length) {
    c(as.numeric(x), rep(NA, times=target_length - length(x)))
  } else {
    as.numeric(x)
  }
}

to_granges <- function(genelist) {
  GRanges(ranges   = IRanges(start=genelist$metagene_start, end=genelist$metagene_stop), 
          seqnames = genelist$chr)
}

get_read_matrix <- function(genelist, sample.cov, pad_length, reverse=FALSE) {
  reads.gr   <- to_granges(genelist)
  reads.rl   <- as(reads.gr, "RangesList")
  reads.view <- RleViewsList(rleList=sample.cov[names(reads.rl)], rangesList=reads.rl)
  reads      <- viewApply(reads.view, function(x) pad_with_nas(x, pad_length, reverse=reverse))
  reads.m    <- matrix(unlist(sapply(reads, sapply, as.numeric)), nrow=nrow(genelist), byrow=TRUE)
  reads.m
}

get_metagene_reads <- function(sample.cov, gene_list, 
                               before_tss=200, after_tss=1500, 
                               smooth=NULL, sample=NULL, 
                               normalization_target=NULL,
                               confidence_intervals=FALSE,
                               return_read_matrix=FALSE) {

  sample_object_name <- deparse(substitute(sample.cov))

  # update the provided object in case it was built with a previous version of IRanges
  sample.cov <- updateObject(sample.cov)

  # total signal in sample
  sample.ts <- calc_total_signal_for_coverage_object(sample.cov)
  
  # normalization factor
  if(is.null(normalization_target)) {
    normalization.factor <- 1
  } else {
    normalization.factor <- normalization_target / sample.ts
  }

  message("Metagene call for ", nrow(gene_list), " genes (", sample_object_name, ")")

  check_gene_list_for_proper_columns(gene_list)

  chr_lengths <- sapply(sample.cov, length)
  
  genes.p <- subset(gene_list, strand == 1)
  genes.n <- subset(gene_list, strand == -1)
  
  reads.p <- NULL
  reads.n <- NULL
  
  if(nrow(genes.p) > 0) {
    genes.p <- transform(genes.p,  metagene_start = start - before_tss,
                                   metagene_stop  = pmin(stop, start + after_tss))
    original_count <- nrow(genes.p)
    genes.p <- subset(genes.p, metagene_start > 0 & metagene_stop <= chr_lengths[as.character(chr)])
    if(nrow(genes.p) < original_count) {
      message("Removing ", original_count - nrow(genes.p), " gene(s) on positive strand due to chromosome boundary")
    }
    reads.p <- get_read_matrix(genes.p, sample.cov, pad_length=before_tss + after_tss + 1, reverse=FALSE)
  }

  if(nrow(genes.n) > 0) {
    genes.n <- transform(genes.n,  metagene_start = pmax(start, stop - after_tss),
                                   metagene_stop  = stop + before_tss)
    original_count <- nrow(genes.n)
    genes.n <- subset(genes.n, metagene_start > 0 & metagene_stop <= chr_lengths[as.character(chr)])
    if(nrow(genes.n) < original_count) {
      message("Removing ", original_count - nrow(genes.n), " gene(s) on negative strand due to chromosome boundary")
    }
    reads.n <- get_read_matrix(genes.n, sample.cov, pad_length=before_tss + after_tss + 1, reverse=TRUE)
  }
  
  reads <- rbind(reads.p, reads.n)
        
  results <- data.frame(tss_distance=(-1*before_tss):after_tss, avg_reads=NA)
  results$avg_reads <- colMeans(reads, na.rm=T) * normalization.factor

  if(confidence_intervals) {
    errors <- qnorm(0.975) * colSds(reads, na.rm=T) / sqrt(nrow(reads))
    results$std_dev <- colSds(reads, na.rm=T)
    results <- transform(results, ci_upper = avg_reads + errors,
                                  ci_lower = avg_reads - errors)
    if(!is.null(smooth)) {
      results$ci_upper_smooth <- as.numeric(runmean(Rle(results$ci_upper), smooth, "constant"))
      results$ci_lower_smooth <- as.numeric(runmean(Rle(results$ci_lower), smooth, "constant"))
    }
  }
  
  if(!is.null(smooth)) results$smooth <- as.numeric(runmean(Rle(results$avg_reads), smooth, "constant"))
  if(!is.null(sample)) results$sample <- sample
  if(return_read_matrix) {
    list(results=results, read_matrix=reads)
  } else {
    results
  }
}

get_metagene_enrichment <- function(ip.cov, bg.cov, gene_list, before_tss=200, after_tss=1500, 
                                    smooth=NULL, sample=NULL, confidence_intervals=FALSE,
                                    return_read_matrix=FALSE) {

  ip_object_name <- deparse(substitute(ip.cov))
  bg_object_name <- deparse(substitute(bg.cov))

  ip.cov <- updateObject(ip.cov)
  bg.cov <- updateObject(bg.cov)

  ip.sum <- calc_total_signal_for_coverage_object(ip.cov)
  bg.sum <- calc_total_signal_for_coverage_object(bg.cov)

  message("Metagene enrichment call for ", nrow(gene_list), " genes (", ip_object_name, " over ", bg_object_name, ")")

  check_gene_list_for_proper_columns(gene_list)

  chr_lengths <- sapply(ip.cov, length)
  
  genes.p <- subset(gene_list, strand == 1)
  genes.n <- subset(gene_list, strand == -1)

  reads.ip.p <- NULL
  reads.bg.p <- NULL
  
  if(nrow(genes.p) > 0) {
    genes.p <- transform(genes.p,  metagene_start = start - before_tss,
                                   metagene_stop  = pmin(stop, start + after_tss))
    original_count <- nrow(genes.p)
    genes.p <- subset(genes.p, metagene_start > 0 & metagene_stop <= chr_lengths[as.character(chr)])
    if(nrow(genes.p) < original_count) {
      message("Removing ", original_count - nrow(genes.p), " gene(s) on positive strand due to chromosome boundary")
    }
    reads.ip.p <- get_read_matrix(genes.p, ip.cov, pad_length=before_tss + after_tss + 1, reverse=FALSE)
    reads.bg.p <- get_read_matrix(genes.p, bg.cov, pad_length=before_tss + after_tss + 1, reverse=FALSE)
  }

  reads.ip.n <- NULL
  reads.bg.n <- NULL

  if(nrow(genes.n) > 0) {
    genes.n <- transform(genes.n,  metagene_start = pmax(start, stop - after_tss),
                                   metagene_stop  = stop + before_tss)
    original_count <- nrow(genes.n)
    genes.n <- subset(genes.n, metagene_start > 0 & metagene_stop <= chr_lengths[as.character(chr)])
    if(nrow(genes.n) < original_count) {
      message("Removing ", original_count - nrow(genes.n), " gene(s) on negative strand due to chromosome boundary")
    }
    reads.ip.n <- get_read_matrix(genes.n, ip.cov, pad_length=before_tss + after_tss + 1, reverse=TRUE)
    reads.bg.n <- get_read_matrix(genes.n, bg.cov, pad_length=before_tss + after_tss + 1, reverse=TRUE)
  }
  
  reads.ip <- rbind(reads.ip.p, reads.ip.n)
  reads.bg <- rbind(reads.bg.p, reads.bg.n)

  reads.e <- (reads.ip / ip.sum) / (reads.bg / bg.sum)
  reads.e[!is.finite(reads.e)] <- NA

  results <- data.frame(tss_distance=(-1*before_tss):after_tss, avg_reads_ip=NA, avg_reads_bg=NA)

  results$avg_reads_ip <- colMeans(reads.ip, na.rm=T) / ip.sum
  results$avg_reads_bg <- colMeans(reads.bg, na.rm=T) / bg.sum

  results$enrichment <- results$avg_reads_ip / results$avg_reads_bg
  #results$enrichment <- colMeans(reads.e, na.rm=T)
  
  if(confidence_intervals) {
    errors <- qnorm(0.975) * colSds(reads.e, na.rm=T) / sqrt(nrow(reads.e))
    results <- transform(results, ci_upper = enrichment + errors,
                                  ci_lower = enrichment - errors)
    if(!is.null(smooth)) results$ci_upper_smooth <- as.numeric(runmean(Rle(results$ci_upper), smooth, "constant"))
    if(!is.null(smooth)) results$ci_lower_smooth <- as.numeric(runmean(Rle(results$ci_lower), smooth, "constant"))
  }

  if(!is.null(smooth)) results$smooth <- as.numeric(runmean(Rle(results$enrichment), smooth, "constant"))
  if(!is.null(sample)) results$sample <- sample

  if(return_read_matrix) {
    list(results=results, read_matrix=reads.e)
  } else {
    results
  }
}


