##############################################
# Plot tandem repeats from SRF
##############################################
load("centromereCandidates.RData")
load("chromSizes.RData")

list.dirs("Annelida_SRF/")
species <- c("NNAJ", "CMAT", "EAND", "LTER", "MVUL", "WPIG", "HNIP", "HMAN", "TLAP", "SBEN", "PROT", "PECH", "AVIR", "APAC", "SNUD", "OFUS", "PMAX", "LLON")

SRFresults <- lapply(species, function(x){
  w <- read.delim(paste("Annelida_SRF/", x, "/srf.bed", sep = ""), header = FALSE)
})
names(SRFresults) <- species

MVULlevels <- rev(c("MVUL.MV040", "MVUL.MV038", "MVUL.MV035", "MVUL.MV033", "MVUL.MV039", "MVUL.MV031", "MVUL.MV037", "MVUL.MV030", "MVUL.MV036", "MVUL.MV028", "MVUL.MV029", "MVUL.MV026", "MVUL.MV027", "MVUL.MV022", "MVUL.MV032", "MVUL.MV021", "MVUL.MV024", "MVUL.MV020", "MVUL.MV023", "MVUL.MV018", "MVUL.MV016", "MVUL.MV012", "MVUL.MV017", "MVUL.MV010", "MVUL.MV015", "MVUL.MV009", "MVUL.MV014", "MVUL.MV008", "MVUL.MV019", "MVUL.MV007", "MVUL.MV011", "MVUL.MV006", "MVUL.MV013", "MVUL.MV004", "MVUL.MV034", "MVUL.MV003", "MVUL.MV002", "MVUL.MV041", "MVUL.MV025", "MVUL.MV005", "MVUL.MV001"))
CMATlevels <- paste("CMAT.chr_", 1:17, sep = "")
NNAJlevels <- paste("NNAJ.chr_", 1:17, sep = "")
LLONlevels <- paste("LLON.", 1:19, sep = "")
PECHlevels <- paste("PECH.Superscaffold", 1:14, sep = "")
SBENlevels <- paste("SBEN.Ch_", 1:11, sep = "")
EANDlevels <- sizes$EAND$names[1:11]
HMANlevels <- sizes$HMAN$names[1:13]
HNIPlevels <- sizes$HNIP$names[1:11]
OFUSlevels <- sizes$OFUS$names[1:12]
SNUDlevels <- sizes$SNUD$names[1:17]
WPIGlevels <- sizes$WPIG$names[1:11]
PMAXlevels <- sizes$PMAX$names[1:19]

toPlotTRs <- function(speciesID){
  x <- which(names(SRFresults) == speciesID)
  temp <- SRFresults[[x]]
  if(speciesID == "SNUD")
    temp[, 1] <- paste("SNUD.", temp[, 1], sep = "")
  temp <- temp[temp[, 1] %in% sizes[[which(names(sizes) == speciesID)]]$names[sizes[[which(names(sizes) == speciesID)]]$size > minChromSize], ]
  if(nrow(temp) > 0){
    TRs <- GRanges(seqnames = temp[, 1], ranges = IRanges(start = temp[, 2], end = temp[, 3]), strand = rep("*", nrow(temp)), family = temp[, 4])
    seqlengths(TRs) <- sapply(unique(seqnames(TRs)), function(y){
      sizes[[which(names(sizes) == speciesID)]]$size[sizes[[which(names(sizes) == speciesID)]]$names == y]
    })
    test <- lapply(unique(seqnames(TRs)), function(z){
      covs <- coverage(TRs)
      tempcov <- covs[names(covs) == z][[1]]
      end <- bwidth * floor(length(tempcov)/bwidth)
      windows <- IRanges(start=seq(1, end, bwidth), width = bwidth)
      cov_by_wnd <- Views(tempcov, windows)
      finalCov <- viewMeans(cov_by_wnd)
      reg <- data.frame(seqnames = z, start = start(windows), end = end(windows), coverage = finalCov * 100)
    })
    cand <- do.call("rbind", test)
  } else{
    cand <- NULL
  }
  return(cand)
}

plotTRsinChroms <- function(species){
  sp1 <- species
  temp <- sizes[[which(names(sizes) == sp1)]]
  temp <- temp[temp$size > minChromSize, ]
  summary <- data.frame(names = factor(temp$names), sizes = temp$size)
  if(sp1 %in% c("MVUL", "CMAT", "NNAJ", "LLON", "PECH", "SBEN", "HMAN", "HNIP", "WPIG", "EAND", "SNUD", "OFUS")){
    summary$names <- factor(summary$names, levels = eval(parse(text = paste(sp1, "levels", sep = "")))[eval(parse(text = paste(sp1, "levels", sep = ""))) %in% summary$names])
  } else{
    summary$names <- factor(summary$names)
  }
  temp2 <- unlist(GRangesList(centromereCandidates[[which(names(centromereCandidates) == sp1)]]))
  names(temp2) <- 1:length(temp2)
  centromeres <- as.data.frame(temp2)
  centromeres <- centromeres[centromeres$seqnames %in% temp$names, ]
  if(sp1 %in% c("MVUL", "CMAT", "NNAJ", "LLON", "PECH", "SBEN", "HMAN", "HNIP", "WPIG", "EAND", "SNUD", "OFUS")){
    centromeres$seqnames <- factor(centromeres$seqnames, levels = eval(parse(text = paste(sp1, "levels", sep = "")))[eval(parse(text = paste(sp1, "levels", sep = ""))) %in% temp$names])
  } else{
    centromeres$seqnames <- factor(centromeres$seqnames)
  }
  TRsofInterest <- toPlotTRs(sp1)
  TRsofInterest <- TRsofInterest[TRsofInterest$coverage > 0, ]
  if(!is.null(TRsofInterest)){
    if(sp1 %in% c("MVUL", "CMAT", "NNAJ", "LLON", "PECH", "SBEN", "HMAN", "HNIP", "WPIG", "EAND", "SNUD", "OFUS")){
      TRsofInterest$seqnames <- factor(TRsofInterest$seqnames, levels = eval(parse(text = paste(sp1, "levels", sep = "")))[eval(parse(text = paste(sp1, "levels", sep = ""))) %in% temp$names])
    } else{
      TRsofInterest$seqnames <- factor(TRsofInterest$seqnames)
    }
    ideogram.boxes <- ggplot(summary, aes(x = names, y = sizes)) + 
      geom_rect(
        aes(ymax = sizes), 
        ymin = 1,
        xmin = as.numeric(summary$names) - 0.4,
        xmax = as.numeric(summary$names) + 0.4, 
        fill = "grey95",
        colour = "dimgrey") + 
      scale_y_continuous(
        limits = c(0, ceiling(max(summary$sizes)/1e6) * 1e6),
        labels = label_number(
          accuracy = 1,
          scale = 1e-6,
          suffix = "Mbp"),
        expand = c(0, 100)) + 
      coord_flip(clip = "off") + 
      theme(axis.text.y = element_text(
        colour = "black",
        size = 8,
        margin = margin(r = 5)),
        axis.text.x = element_text(size = 8),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())
    ideogram.TEs <- ideogram.boxes +
      geom_segment(data = TRsofInterest,
                   aes(y = start,
                       yend = start,
                       colour = coverage),
                   x = as.numeric(TRsofInterest$seqnames) - 0.4,
                   xend = as.numeric(TRsofInterest$seqnames) + 0.4,
                   size = 0.7) +
      scale_colour_gradient(
        name=paste("Percentage of bases covered by TRs", sep = ""), 
        limits = c(min(na.omit(TRsofInterest$coverage)), max(na.omit(TRsofInterest$coverage))),
        low = "tan1", high = "mediumorchid4", na.value = "transparent"
      )
    ideogram.centromeres <- ideogram.TEs +
      geom_segment(data = centromeres,
                   aes(y = start,
                       yend = end),
                   x = as.numeric(centromeres$seqnames),
                   xend = as.numeric(centromeres$seqnames),
                   size = 2) +
      theme(legend.position="bottom", legend.box = "horizontal")
    ggsave(paste(sp1, "_TRs.svg", sep = ""), ideogram.centromeres, width = 30, height = 30, units = "cm", bg = "white")
    ggsave(paste(sp1, "_TRs.png", sep = ""), ideogram.centromeres, device = "png", width = 30, height = 30, units = "cm", bg = "white")
  }
}

sapply(species, function(x){
  plotTRsinChroms(x)
})