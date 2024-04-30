##################################################################
#
# To identify centromeres
#
##################################################################

library(ggplot2) 
library(ShortRead)
library(pbmcapply)
library(GenomicRanges)

# Load pre-existing data objects containing chromosome sizes and TE annotation
load("chromSizes.RData")
load("Superfamilies_RModeler.RData")

species <- c("NNAJ", "CMAT", "EAND", "LTER", "MVUL", "WPIG", "HNIP", "HMAN", "TLAP", "SBEN", "PROT", "PECH", "AVIR", "APAC", "SNUD", "OFUS", "PMAX", "LLON")

# Minimum chromosome size threshold in base pairs
minChromSize <- 5000000

# Threshold value for identifying candidate regions, will consider the % of the genome with the highest TE density
thold <- .91

# Bin width for coverage calculations in base pairs
bwidth <- 50000L

# Number of bins on either side of the current bin used for smoothing coverage values
sidebins <- 25

centromereCandidates <- lapply(species, function(x){
  tempTEs <- superfamilies[[which(names(superfamilies) == x)]]
  tempSizes <- sizes[[which(names(sizes) == x)]]
  tempTEs <- tempTEs[tempTEs$qname %in% tempSizes$names[tempSizes$size > minChromSize] & tempTEs$p_sub >= 0, ]
  
  # Create GRanges object representing filtered transposable elements with annotations
  TEs <- GRanges(seqnames = tempTEs$qname, ranges = IRanges(start = tempTEs$qstart, end = tempTEs$qend), strand = rep("*", nrow(tempTEs)), family = tempTEs$tname, superfamily = tempTEs$tclass, order = sapply(strsplit(tempTEs$tclass, "/"), "[[", 1))
  
  # Set sequence lengths for TEs based on chromosome sizes
  seqlengths(TEs) <- sapply(unique(seqnames(TEs)), function(y){
    tempSizes$size[tempSizes$names == y]
  })
  
  # Calculate coverage for transposable elements
  covs <- coverage(TEs)
  
  # Filter chromosomes larger than the minimum size for the current species
  tempsubset <- sizes[[which(names(sizes) == x)]]
  tempsubset <- tempsubset$names[tempsubset$size > minChromSize]
  
  candidateRegions <- lapply(tempsubset, function(y){
    # Get coverage data and calculate end position for binning on the current chromosome
    tempcov <- covs[names(covs) == y][[1]]
    end <- bwidth * floor(length(tempcov)/bwidth)
    
    # Create windows of size bwidth along the chromosome
    windows <- IRanges(start=seq(1, end, bwidth), width = bwidth)
    
    # Calculate average coverage per window and perform smoothing
    cov_by_wnd <- Views(tempcov, windows)
    finalCov <- viewMeans(cov_by_wnd)
    smooth <- sapply(1:length(finalCov), function(z){
      start <- z - sidebins
      if(start < 1) start <- 1
      end <- z + sidebins
      if(end > length(finalCov)) end <- length(finalCov)
      median(finalCov[start:end])
    })
    
    # Identify candidate regions based on smoothed coverage exceeding a threshold
    th <- quantile(smooth, thold)
    reg <- reduce(GRanges(seqnames = rep(y, sum(smooth > th)), windows[smooth > th], strand = rep("*", sum(smooth > th))))
    reg <- reg[width(reg) > 20 * bwidth]
    
    # Set names of candidate regions based on chromosome names
  })
  names(candidateRegions) <- tempsubset
  return(candidateRegions)
})
names(centromereCandidates) <- species
#save(centromereCandidates, file = "centromereCandidates.RData")

##################################################################
#
# Identify TEs enriched in the centromeres
#
##################################################################

TEsInCentromeres <- lapply(species, function(x){
  tempTEs <- superfamilies[[which(names(superfamilies) == x)]]
  tempSizes <- sizes[[which(names(sizes) == x)]]
  tempTEs <- tempTEs[tempTEs$qname %in% tempSizes$names[tempSizes$size > minChromSize] & tempTEs$p_sub >= 0, ]
  TEs <- GRanges(seqnames = tempTEs$qname, ranges = IRanges(start = tempTEs$qstart, end = tempTEs$qend), strand = rep("*", nrow(tempTEs)), family = tempTEs$tname, superfamily = tempTEs$tclass, order = sapply(strsplit(tempTEs$tclass, "/"), "[[", 1))
  seqlengths(TEs) <- sapply(unique(seqnames(TEs)), function(y){
    tempSizes$size[tempSizes$names == y]
  })
  candidateRegions <- centromereCandidates[[which(names(centromereCandidates) == x)]]
  pvalsTEsinbins <- pbmcapply::pbmclapply(1:length(unique(TEs$superfamily)), function(y){
    subsetTEs <- which(TEs$superfamily == unique(TEs$superfamily)[y])
    pvals <- sapply(1:length(candidateRegions), function(z){
      toCount <- findOverlaps(candidateRegions[[z]], TEs)
      a <- sum(subjectHits(toCount) %in% subsetTEs)
      b <- sum(seqnames(TEs) %in% names(candidateRegions)[z] & TEs$superfamily %in% unique(TEs$superfamily)[y]) - a
      c <- length(unique(subjectHits(toCount))) - a
      d <- sum(seqnames(TEs) %in% names(candidateRegions)[z]) - a - b - c
      test <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")$p.value
    })
    adjpvals <- p.adjust(pvals, "BH")
    return(adjpvals)
  }, mc.cores = 10)
  pvalsTEsinbins <- do.call("rbind", pvalsTEsinbins)
  rownames(pvalsTEsinbins) <- unique(TEs$superfamily)
  colnames(pvalsTEsinbins) <- names(candidateRegions)
  superfamily <- sapply(rownames(pvalsTEsinbins), function(y){
    TEs$superfamily[which(TEs$superfamily == y)[1]]
  })
  enrichedSuperfamilies <- table(superfamily[which(rowSums(pvalsTEsinbins < .05) >= trunc(sum(colSums(pvalsTEsinbins)/nrow(pvalsTEsinbins) != 1) * .85))])
  out <- list(enrichedSuperfamilies, pvalsTEsinbins)
  return(out)
})
names(TEsInCentromeres) <- species
