
fisher_test_2x2 <- function(groupA, groupB, gene.universe, verbose=FALSE) {
  if(verbose) message("fisher_test_2x2()")
  groupA            <- unique(as.character(groupA))
  groupB            <- unique(as.character(groupB))
  gene.universe     <- unique(as.character(gene.universe))

  # First, restrict both groups to the provided universe of genes
  groupA <- groupA[groupA %in% gene.universe]
  groupB <- groupB[groupB %in% gene.universe]

  if(length(groupA) == 0 | length(groupB) == 0) {
    message("Zero elements after universe filtering. groupA: ", length(groupA), " groupB: ", length(groupB))
    result.df <- data.frame(test_type="Enrichment", 
                            pvalue=1,
                            overlap = 0,
                            totalA  = length(groupA),
                            totalB  = length(groupB),
                            universe = length(gene.universe),
                            enrichment = 1)
    return(result.df)
  }
  
  # Calculate the 4 values for the contingency table
  count.Aonly   <- length(groupA[!(groupA %in% groupB)])
  count.overlap <- length(groupA[groupA %in% groupB])
  count.neither <- length(gene.universe[!(gene.universe %in% groupA) & !(gene.universe %in% groupB)])
  count.Bonly   <- length(groupB[!(groupB %in% groupA)])
  
  if(verbose) {
    message(" Universe size: ", length(gene.universe))
    message(" A only: ", count.Aonly)
    message(" Overlap A & B: ", count.overlap)
    message(" Neither A nor B: ", count.neither)
    message(" B only: ", count.Bonly)
  }
  
  # Build the 2x2 matrix
  m <- matrix(c(count.overlap, count.Aonly, count.Bonly, count.neither), nrow=2, byrow=T)
  
  # Decide on enrichment test or depletion test
  prop1 <- count.overlap / length(groupA)
  prop2 <- count.Bonly / (length(gene.universe) - length(groupA))
  if(is.na(prop2)) prop2 <- -1
  
  if(verbose) {
    message("Proportion 1: ", ifelse(is.numeric(prop1), round(prop1, 4), "<not numeric>"))
    message("Proportion 2: ", ifelse(is.numeric(prop2), round(prop2, 4), "<not numeric>"))
  }
  
  alt.test <- ifelse(prop1 > prop2, "greater", "less")
  
  # Perform the test
  result <- fisher.test(m, alternative=alt.test)
  
  obs <- count.overlap / length(groupA)
  exp <- length(groupB) / length(gene.universe)
  
  enrichment <- obs / exp
  
  # Return the result
  result.df <- data.frame(test_type=ifelse(alt.test == "greater", "Enrichment", "Depletion"), 
                          pvalue=result$p.value,
                          overlap = as.integer(count.overlap),
                          totalA  = length(groupA),
                          totalB  = length(groupB),
                          universe = length(gene.universe),
                          enrichment = enrichment)
  result.df
}

