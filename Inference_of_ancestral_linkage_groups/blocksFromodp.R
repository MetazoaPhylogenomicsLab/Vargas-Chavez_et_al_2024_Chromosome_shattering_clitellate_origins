load("Superfamilies_RModeler.RData")
load("centromereCandidates.RData")
load("chromSizes.RData")
species <- c("CMAT", "EAND", "HMAN", "HNIP", "LLON", "MVUL", "NNAJ", "OFUS", "PECH", "PMAX", "SBEN", "SNUD", "WPIG")

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

getBlocks <- function(speciesID, directory, threshold, LGsID){
  sp1 <- speciesID
  print(sp1)
  temp <- read.delim(directory)
  temp <- temp[temp[, 8] %in% eval(parse(text = paste(sp1, "levels", sep = ""))), ]
  temp <- temp[!duplicated(temp[, 7]), ]
  temp <- temp[temp[, 16] < 0.05, ]
  toPlot <- NA
  if(length(temp[, 8]) > 0 & all(table(temp[, 8]) > 10)){
    colorsTemp <- temp[!duplicated(temp$color), 2:3]
    colors <- unlist(colorsTemp$color)
    names(colors) <- colorsTemp$gene_group
    colors <- colors[order(names(colors))]
    rownames(temp) <- temp[, 7]
    temp <- temp[, c(9, 2, 8)]
    colnames(temp) <- c("pos", "match", "chrom")
    gff <- read.delim(paste("/mnt/disk2/Earthworms/References/processed/", sp1, ".gff", sep = ""), header = FALSE, sep = "\t")
    gff[, 9] <- sub(":", "_", gff[, 9])
    gff <- gff[gff[, 3] == "gene", ]
    tempS <- sizes[[which(names(superfamilies) == sp1)]]
    chroms2plot <- unique(temp$chrom)[unique(temp$chrom) %in% tempS$names[tempS$size > 5000000]]
    toPlot <- lapply(chroms2plot, function(y){
      print(y)
      temp2 <- temp[temp$chrom == y, ]
      tempgff <- gff[gff[, 1] == y, ]
      tempgff <- tempgff[order(tempgff[, 4]), ]
      temp2$pos <- sapply(rownames(temp2), function(z){
        which(tempgff[, 9] == z)
      })
      temp2 <- temp2[order(temp2$pos), ]
      pos <- Rle(temp2$match)
      cands <- which(runLength(pos) <= threshold)
      cands <- cands[-c(1, length(cands))]
      remove <- cands[which(runValue(pos)[cands - 1] == runValue(pos)[cands + 1])]
      if(length(remove) > 0){
        cands <- cands[-remove]
        toRemove <- unlist(lapply(cands, function(z){
          out <- (cumsum(runLength(pos))[z] - runLength(pos)[z] + 1):cumsum(runLength(pos))[z]
        }))
        temp2 <- temp2[-toRemove, ]
      }
      pos <- rle(temp2$match)
      pos2 <- cumsum(pos$lengths)
      starts <- rownames(temp2)[c(0, pos2[1:length(pos2) - 1]) + 1]
      ends <- rownames(temp2)[pos2]
      byChrom <- lapply(1:length(pos2), function(z){
        data.frame(start = gff[which(gff[, 9] == starts[z]), 4], end = gff[which(gff[, 9] == ends[z]), 5], match = pos$values[z], chrom = y)
      })
      byChrom <- do.call("rbind", byChrom)
    })
    toPlot <- do.call("rbind", toPlot)
    if(sp1 %in% c("MVUL", "CMAT", "NNAJ", "LLON", "PECH", "SBEN", "HMAN", "HNIP", "WPIG", "EAND", "SNUD", "OFUS")){
      toPlot$chrom <- factor(toPlot$chrom, levels = eval(parse(text = paste(sp1, "levels", sep = "")))[eval(parse(text = paste(sp1, "levels", sep = ""))) %in% toPlot$chrom])
    } else{
      toPlot$chrom <- factor(toPlot$chrom)
    }
    write.table(toPlot[, c(6, 1, 2, 5)], paste(sp1, "_", LGsID, ".txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  return(toPlot)
}

blocks <- lapply(species, function(x){
  leech <- getBlocks(x, paste("odp/step2-figures/ALG-species_plots/odp_LGs_from_HMAN_HNIP_WPIG_", x, "_xy_reciprocal_best_hits.plotted.rbh", sep = ""), 1, "LeechLGs")
  list(leech)
})
