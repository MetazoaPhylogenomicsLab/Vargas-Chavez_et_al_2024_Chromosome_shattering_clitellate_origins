##########################
# Ribbon plots
##########################
library(BSgenome)
library(dplyr)
library(igraph)
library(ggplot2)
library(ggnewscale)
library(ggpubr)
library(parallel)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(RIdeogram)
library(scales)
library(ShortRead)
library(tgutil)
library(tibble)

withAnnot <- c("CMAT", "EAND", "ECRY", "HMAN", "HNIP", "LLON", "MVUL", "NNAJ", "OFUS", "PECH", "PMAX", "SBEN", "SNUD", "WPIG")
load("chromSizes.RData")
minChromSize <- 5000000
bwidth <- 50000L
th <- 0.91

karyotypes <- lapply(which(names(sizes) %in% withAnnot), function(x){
  karyotype <- data.frame(Chr = sizes[[x]][, 1], Start = 1, End = sizes[[x]][, 2], fill = "ffffffff", species = names(sizes)[x], size = 12, color = "252525")
  karyotype <- karyotype[karyotype$End > minChromSize, ]
  karyotype <- karyotype[order(karyotype$End, decreasing = TRUE), ]
  karyotype_mapping <- setNames(1:nrow(karyotype), karyotype$Chr)
  karyotype$orig_ids <- karyotype$Chr
  karyotype$Chr <- karyotype_mapping[karyotype$Chr]
  return(karyotype)
})
names(karyotypes) <- names(sizes)[which(names(sizes) %in% withAnnot)]

getOrder <- function(sp1, sp2, newOrder = NULL, path = "odp/step2-figures/synteny_coloredby_BCnS_LGs_merged_Trochozoa"){
  temprbh <- read.delim(paste(path, "/", sort(c(sp1, sp2))[1], "_", sort(c(sp1, sp2))[2], "_xy_reciprocal_best_hits.", unlist(strsplit(path, "/synteny_"))[2], ".plotted.rbh", sep = ""))
  sp1karyotype <- karyotypes[[which(names(karyotypes) == sp1)]]
  if(!is.null(newOrder))
    sp1karyotype <- sp1karyotype[newOrder, ]
  sp2karyotype <- karyotypes[[which(names(karyotypes) == sp2)]]
  temprbh <- temprbh[get(paste(sp1, "_scaf", sep = ""), temprbh) %in% sp1karyotype$orig_ids & get(paste(sp2, "_scaf", sep = ""), temprbh) %in% sp2karyotype$orig_ids, ]
  usedsp2karyotype <- rep(FALSE, nrow(sp2karyotype))
  sp2Order <- unlist(sapply(1:nrow(sp1karyotype), function(x){
    allMatches <- table(get(paste(sp2, "_scaf", sep = ""), temprbh)[get(paste(sp1, "_scaf", sep = ""), temprbh) == sp1karyotype$orig_ids[x]])
    allMatches <- sort(allMatches, decreasing = TRUE)
    allMatches <- allMatches[allMatches > 15]
    temppos <- sapply(names(allMatches), function(y){
      which(sp2karyotype$orig_ids == y)
    })
    temppos <- temppos[usedsp2karyotype[temppos] == FALSE]
    usedsp2karyotype[temppos] <<- TRUE
    return(temppos)
  }))
  return(sp2Order)
}

getSyntenyPlot <- function(sp1, sp2, order1 = NULL, order2 = NULL, splitbysp1Chrom = FALSE, path = "odp/step2-figures/synteny_coloredby_BCnS_LGs_merged_Trochozoa", flip1 = NULL, flip2 = NULL){
  kar1 <- karyotypes[[which(names(karyotypes) == sp1)]]
  if(!is.null(order1))
    kar1 <- kar1[order1, ]
  kar1_remap <- kar1$Chr
  names(kar1_remap) <- kar1$orig_ids
  kar2 <- karyotypes[[which(names(karyotypes) == sp2)]]
  if(!is.null(order2))
    kar2 <- kar2[order2, ]
  kar2_remap <- kar2$Chr
  names(kar2_remap) <- kar2$orig_ids
  karyotype_all <- rbind(kar1, kar2, make.row.names = FALSE)
  karyotype_all$orig_ids <- NULL
  karyotype_all <- as.data.frame(karyotype_all, stringsAsFactors = FALSE)
  temprbh <- read.delim(paste(path, "/", sort(c(sp1, sp2))[1], "_", sort(c(sp1, sp2))[2], "_xy_reciprocal_best_hits.", unlist(strsplit(path, "/synteny_"))[2], ".plotted.rbh", sep = ""))
  synteny_all <- temprbh[get(paste(sp1, "_scaf", sep = ""), temprbh) %in% kar1$orig_ids & get(paste(sp2, "_scaf", sep = ""), temprbh) %in% kar2$orig_ids & temprbh$whole_FET < 0.05, ]
  if(nrow(synteny_all) < 200){
    temprbh <- read.delim(paste(path, "/", sort(c(sp1, sp2))[1], "_", sort(c(sp1, sp2))[2], "_xy_reciprocal_best_hits.", unlist(strsplit(path, "/synteny_"))[2], ".plotted.rbh", sep = ""))
    synteny_all <- temprbh[get(paste(sp1, "_scaf", sep = ""), temprbh) %in% kar1$orig_ids & get(paste(sp2, "_scaf", sep = ""), temprbh) %in% kar2$orig_ids, ]
  }
  synteny_all <- data.frame(Species_1 = get(paste(sp1, "_scaf", sep = ""), synteny_all), Start_1 = get(paste(sp1, "_pos", sep = ""), synteny_all), End_1 = get(paste(sp1, "_pos", sep = ""), synteny_all) + 1000, Species_2 = get(paste(sp2, "_scaf", sep = ""), synteny_all), Start_2 = get(paste(sp2, "_pos", sep = ""), synteny_all), End_2 = get(paste(sp2, "_pos", sep = ""), synteny_all) + 1000, fill = sub("#", "", synteny_all$color))
  synteny_all$Species_1 <- kar1_remap[match(synteny_all$Species_1, names(kar1_remap))]
  synteny_all$Species_2 <- kar2_remap[match(synteny_all$Species_2, names(kar2_remap))]
  synteny_all <- synteny_all %>%
    group_by(Species_1, Species_2) %>%
    filter(n() > 5) %>%
    ungroup()
  synteny_all <- synteny_all[with(synteny_all, order(fill)),]
  synteny_all <- as.data.frame(synteny_all)
  if(!is.null(order1))
    synteny_all$Species_1 <- match(synteny_all$Species_1, kar1_remap)
  if(!is.null(order2))
    synteny_all$Species_2 <- match(synteny_all$Species_2, kar2_remap)
  synteny_all$fill <- sub("000000", "cccccc", synteny_all$fill)
  synteny_all <- synteny_all[c(which(synteny_all$fill == "cccccc"), sample(which(synteny_all$fill != "cccccc"), sum(synteny_all$fill != "cccccc"))), ]
  synteny_filt <- synteny_all[synteny_all$fill != "cccccc", ]
  synteny_filt <- synteny_filt %>%
    group_by(Species_1, Species_2) %>%
    filter(n() > 5) %>%
    ungroup()
  synteny_filt <- synteny_filt[with(synteny_filt, order(fill)),]
  synteny_filt <- as.data.frame(synteny_filt)
  synteny_filt <- synteny_filt[sample(which(synteny_filt$fill != "cccccc"), sum(synteny_filt$fill != "cccccc")), ]
  if(!is.null(flip1)){
    sapply(flip1, function(x){
      size <- karyotype_all$End[which(karyotype_all$species == unique(karyotype_all$species)[1])[x]]
      synteny_all$Start_1[synteny_all$Species_1 == x] <<- abs(synteny_all$Start_1[synteny_all$Species_1 == x] - size - 1)
      synteny_all$End_1[synteny_all$Species_1 == x] <<- abs(synteny_all$End_1[synteny_all$Species_1 == x] - size - 1)
      synteny_filt$Start_1[synteny_filt$Species_1 == x] <<- abs(synteny_filt$Start_1[synteny_filt$Species_1 == x] - size - 1)
      synteny_filt$End_1[synteny_filt$Species_1 == x] <<- abs(synteny_filt$End_1[synteny_filt$Species_1 == x] - size - 1)
    })
  }
  if(!is.null(flip2)){
    sapply(flip2, function(x){
      size <- karyotype_all$End[which(karyotype_all$species == unique(karyotype_all$species)[2])[x]]
      synteny_all$Start_2[synteny_all$Species_2 == x] <<- abs(synteny_all$Start_2[synteny_all$Species_2 == x] - size - 1)
      synteny_all$End_2[synteny_all$Species_2 == x] <<- abs(synteny_all$End_2[synteny_all$Species_2 == x] - size - 1)
      synteny_filt$Start_2[synteny_filt$Species_2 == x] <<- abs(synteny_filt$Start_2[synteny_filt$Species_2 == x] - size - 1)
      synteny_filt$End_2[synteny_filt$Species_2 == x] <<- abs(synteny_filt$End_2[synteny_filt$Species_2 == x] - size - 1)
    })
  }
  if(splitbysp1Chrom == FALSE){
    ideogram(karyotype = karyotype_all, synteny = synteny_all, output = paste0("ideograms/svg/", sp1, "_", sp2, "_", unlist(strsplit(path, "/synteny_"))[2], "_synteny.svg"))
    convertSVG(paste0("ideograms/svg/", sp1, "_", sp2, "_", unlist(strsplit(path, "/synteny_"))[2], "_synteny.svg"), file = paste0("ideograms/pdf/", sp1, "_", sp2, "_", unlist(strsplit(path, "/synteny_"))[2], "_synteny.pdf"), device = "pdf")
    ideogram(karyotype = karyotype_all, synteny = synteny_filt, output = paste0("ideograms/svg/", sp1, "_", sp2, "_", unlist(strsplit(path, "/synteny_"))[2], "_synteny_onlyLGs.svg"))
    convertSVG(paste0("ideograms/svg/", sp1, "_", sp2, "_", unlist(strsplit(path, "/synteny_"))[2], "_synteny_onlyLGs.svg"), file = paste0("ideograms/pdf/", sp1, "_", sp2, "_", unlist(strsplit(path, "/synteny_"))[2], "_synteny_onlyLGs.pdf"), device = "pdf")
  } else{
    sapply(karyotype_all$Chr[karyotype_all$species == sp1], function(x){
      ideogram(karyotype = karyotype_all, synteny = synteny_all[synteny_all$Species_1 == x, ], output = paste0("ideograms/svg/", sp1, "_", sp2, "_", unlist(strsplit(path, "/synteny_"))[2], "_synteny_", x, ".svg"))
      convertSVG(paste0("ideograms/svg/", sp1, "_", sp2, "_", unlist(strsplit(path, "/synteny_"))[2], "_synteny_", x, ".svg"), file = paste0("ideograms/pdf/", sp1, "_", sp2, "_", unlist(strsplit(path, "/synteny_"))[2], "_synteny_", x, ".pdf"), device = "pdf")
      if(nrow(synteny_filt[synteny_filt$Species_1 == x, ]) > 0){
        ideogram(karyotype = karyotype_all, synteny = synteny_filt[synteny_filt$Species_1 == x, ], output = paste0("ideograms/svg/", sp1, "_", sp2, "_", unlist(strsplit(path, "/synteny_"))[2], "_synteny_", x, "_onlyLGs.svg"))
        convertSVG(paste0("ideograms/svg/", sp1, "_", sp2, "_", unlist(strsplit(path, "/synteny_"))[2], "_synteny_", x, "_onlyLGs.svg"), file = paste0("ideograms/pdf/", sp1, "_", sp2, "_", unlist(strsplit(path, "/synteny_"))[2], "_synteny_", x, "_onlyLGs.pdf"), device = "pdf")
      }
    })
  }
}

orderA <- getOrder("LLON", "PMAX")
getSyntenyPlot("LLON", "PMAX", order1 = NULL, order2 = orderA)
orderB <- getOrder("PMAX", "OFUS", orderA)
getSyntenyPlot("PMAX", "OFUS", order1 = orderA, order2 = orderB)
orderC <- getOrder("OFUS", "SNUD", orderB)
getSyntenyPlot("OFUS", "SNUD", order1 = orderB, order2 = orderC)
orderD <- getOrder("SNUD", "PECH", orderC)
getSyntenyPlot("SNUD", "PECH", order1 = orderC, order2 = orderD)
orderE <- getOrder("PECH", "SBEN", orderD)
getSyntenyPlot("PECH", "SBEN", order1 = orderD, order2 = orderE)
orderF <- getOrder("SBEN", "HMAN", orderE)
getSyntenyPlot("SBEN", "HMAN", order1 = orderE, order2 = orderF)
orderG <- getOrder("HMAN", "HNIP", orderF)
getSyntenyPlot("HMAN", "HNIP", order1 = orderF, order2 = orderG, flip2 = c(3, 5, 8, 9, 10))
orderH <- getOrder("HNIP", "WPIG", orderG)
getSyntenyPlot("HNIP", "WPIG", order1 = orderG, order2 = orderH, flip1 = c(3, 5, 8, 9, 10), flip2 = c(1, 4, 6, 7, 9, 10, 11))
orderI <- getOrder("WPIG", "NNAJ", orderH)
getSyntenyPlot("WPIG", "NNAJ", order1 = orderH, order2 = orderI, flip1 = c(1, 4, 6, 7, 9, 10, 11))
orderJ <- getOrder("NNAJ", "CMAT", orderI)
getSyntenyPlot("NNAJ", "CMAT", order1 = orderI, order2 = orderJ[c(1:9, 11, 12, 10, 13:17)], flip2 = c(1, 2, 4, 5, 7, 8, 11, 13, 15))
orderK <- getOrder("CMAT", "EAND", orderJ)
getSyntenyPlot("CMAT", "EAND", order1 = orderJ[c(1:9, 11, 12, 10, 13:17)], order2 = orderK, flip1 = c(1, 2, 4, 5, 7, 8, 11, 13, 15), flip2 = c(1, 4, 6, 10, 11))
orderL <- getOrder("EAND", "MVUL", orderK)
getSyntenyPlot("EAND", "MVUL", order1 = orderK, order2 = orderL, flip1 = c(1, 4, 6, 10, 11))

##################
# Oxford dotplot
##################

PMAXvsSNUD <- read.delim("odp/step2-figures/synteny_coloredby_BCnS_LGs_merged_Trochozoa/PMAX_SNUD_xy_reciprocal_best_hits.coloredby_BCnS_LGs_merged_Trochozoa.plotted.rbh")
colors <- PMAXvsSNUD[, c(4, 23)]
colors <- colors[!duplicated(colors[, 1]), ]
colors <- colors[order(colors[, 1]), ]
colors <- colors[c(1:18, 20:21, 19), ]
colors[, 2] <- toupper(colors[, 2])
colorsALGs <- colors[, 2]
names(colorsALGs) <- colors[, 1]
colorsALGs <- sub("#000000", "#cccccc", colorsALGs)

getDotPlot <- function(sp1, sp2, path = "odp/step2-figures/synteny_coloredby_BCnS_LGs_merged_Trochozoa"){
  rbh <- read.delim(paste(path, "/", sort(c(sp1, sp2))[1], "_", sort(c(sp1, sp2))[2], "_xy_reciprocal_best_hits.", unlist(strsplit(path, "/synteny_"))[2], ".plotted.rbh", sep = ""))
  rbh <- rbh[get(paste(sp1, "_scaf", sep = ""), rbh) %in% sizes[[which(names(sizes) == sp1)]][which(sizes[[which(names(sizes) == sp1)]][, 2] > 2000000), 1] & get(paste(sp2, "_scaf", sep = ""), rbh) %in% sizes[[which(names(sizes) == sp2)]][which(sizes[[which(names(sizes) == sp2)]][, 2] > 2000000), 1], ]
  if(sp1 == "PMAX"){
    PMAXorder <- c("PMAX.NC_047026.1", "PMAX.NC_047021.1", "PMAX.NC_047031.1", "PMAX.NC_047018.1", "PMAX.NC_047017.1", "PMAX.NC_047015.1", "PMAX.NC_047032.1", "PMAX.NC_047016.1", "PMAX.NC_047025.1", "PMAX.NC_047022.1", "PMAX.NC_047019.1", "PMAX.NC_047024.1", "PMAX.NC_047023.1", "PMAX.NC_047020.1", "PMAX.NC_047029.1", "PMAX.NC_047028.1", "PMAX.NC_047030.1", "PMAX.NC_047033.1", "PMAX.NC_047027.1")
    rbh <- rbh[order(rbh$PMAX_pos, decreasing = FALSE), ]
    newOrder <- unlist(lapply(PMAXorder, function(x){
      which(rbh$PMAX_scaf == x)
    }))
    rbh <- rbh[newOrder, ]
    rbh$PMAX_plotindex <- 0:(nrow(rbh)- 1)
  }
  rbh <- rbh[c(which(rbh$gene_group == "None"), which(rbh$gene_group != "None")), ]
  sp1chroms <- c(0, sapply(unique(get(paste(sp1, "_scaf", sep = ""), rbh)), function(x){
    max(get(paste(sp1, "_plotindex", sep = ""), rbh)[which(get(paste(sp1, "_scaf", sep = ""), rbh) == x)])
  }))
  sp2chroms <- c(0, sapply(unique(get(paste(sp2, "_scaf", sep = ""), rbh)), function(x){
    max(get(paste(sp2, "_plotindex", sep = ""), rbh)[which(get(paste(sp2, "_scaf", sep = ""), rbh) == x)])
  }))
  rbh$gene_group <- factor(rbh$gene_group, levels = colors[, 1])
  rbh$sp1 <- get(paste(sp1, "_plotindex", sep = ""), rbh)
  rbh$sp2 <- get(paste(sp2, "_plotindex", sep = ""), rbh)
  p1 <- ggplot(rbh, aes(x = sp1, y = sp2, color = gene_group)) + 
    geom_hline(yintercept = sp2chroms, color = "gray", linewidth = 0.5) + 
    geom_vline(xintercept = sp1chroms, color = "gray", linewidth = 0.5) + 
    geom_point(shape = 20) +
    guides(color = guide_legend(title = "BCnS ancestral linkage groups")) +
    scale_color_manual(values = colorsALGs) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(expand = c(0, 0), breaks = NULL, labels = NULL) + 
    scale_y_continuous(expand = c(0, 0), breaks = NULL, labels = NULL) +
    labs(x = sp1, y = sp2) + 
    theme_light()
  return(p1)
}

p1 <- getDotPlot("PMAX", "SNUD")
p2 <- getDotPlot("PMAX", "NNAJ")
p3 <- getDotPlot("PMAX", "ECRY")
p4 <- getDotPlot("HMAN", "SNUD")
p5 <- getDotPlot("HMAN", "NNAJ")
p6 <- getDotPlot("HMAN", "ECRY")

dotplots <- ggarrange(plotlist = list(p1, p4, p2, p5, p3, p6), nrow = 3, ncol = 2, common.legend = TRUE, legend = "right")


