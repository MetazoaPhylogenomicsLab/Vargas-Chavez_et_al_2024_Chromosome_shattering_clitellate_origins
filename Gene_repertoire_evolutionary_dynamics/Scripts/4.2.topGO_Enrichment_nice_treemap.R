### Running functional enrichment using TopGO
# https://bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO.pdf

#Script wrote by Lisandra Benítez Álvarez and Gemma I. Martínez-Redondo on February 2024

#Installing libraries:

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   + install.packages("BiocManager")
# BiocManager::install("topGO")

#install.packages("scico")
#install.packages("reticulate")
#install.packages('rJava')
#install.packages("rjson")
#install.packages("PATH_to/rJython_0.0-4.tar.gz", repos = NULL, type = "source") #https://cran.r-project.org/src/contrib/Archive/rJython/ Archived package
#devtools::install_github("FedericoComoglio/rSymPy")
#BiocManager::install("simplifyEnrichment")

# library(topGO)
# library(dplyr)
# library(data.table)
# library(stringr)
# library(httr)
# library(ggplot2)
# library(scico)
# library(treemap)
# library(reticulate)
# library(rJava)
# library(rjson)
# library(rJython)
# library(rSymPy)
# sympyStart()
# library(GOfuncR)
# library(igraph)
# library(simplifyEnrichment)
# library(ggtree)
# library(DBI)
# library(GOfuncR)
# library(clusterProfiler)

#Set group for analysis
GROUPS <- c("node58", "node50",  "node59", "node60", "node61", "node62", "node63", "node64", "node65", "node66", "node67", "node68", "node69", "node70", "node71", "node72") #list of folders

for (GROUP in GROUPS) {

  sets <- c("LOST", "GAINED", "RETAINED", "DUPLICATED", "EXPANDED") #list of sets
  
  for (set in sets) {
    #set GO categories list to analyse: bpo, mfo, or cco
    categories <- c("BP", "MF", "CC")
    
    #Start loop in categories for set
    for (cat in categories) {
      
      ont = cat
      
      #Set Global directory (GD) and Work directory (WD)
      GD <- paste0("PATH_to_Dir/", GROUP)
      anals <- paste0("Enrichm_", GROUP, "_", set, "_", cat)
      dir.create(paste0(GD, "/", anals))
      WD <- paste0(GD, "/", anals) #create your working directory
      setwd(WD)
      
      #Import background
      geneID2GO <- readMappings(file=paste0("Path_to/Background_seqs_", cat, "_topgo.txt"))
      geneNames <- names(geneID2GO)
      head(geneNames)
      
      #Import GOs list to have terms
      annotation = read.csv(file=paste0("Path_to/Background_seqs_", cat, "_list.txt"), sep= '\t', header=FALSE)
      GOs <- annotation$V1
      GO2term <- go2term(GOs)
      
      #Import gene list
      geneIds = readLines(paste0(GD,"/", set, "_", GROUP, "_seqs.txt")) 
      myInterestingGenes <- as.character(geneIds)
      geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
      names(geneList) <- geneNames
      str(geneList)
      
      GOdata <- new("topGOdata", description ="Simple session", ontology = cat, allGenes = geneList, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = geneID2GO)
      GOdata
      
      resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
      pvalues<-score(resultFisher)
      write.table(pvalues, file=paste0("topGO_enrichment_", anals, "_all_pvalues.txt") ,col.names=FALSE,row.names=TRUE,quote=FALSE, sep = "\t")
      
      data <- cbind(names(pvalues),as.data.table(pvalues))
      data_sig <- data[data$pvalues < 0.05]
      countp <- nrow(data_sig)
      write.table(data_sig,file=paste0("topGO_enrichment_", anals, "_", countp, "_sign.txt"), col.names=FALSE,row.names=FALSE,quote=FALSE, sep = "\t")
      data_sig_info <- merge(data_sig, GO2term, by.x = "V1", by.y = "go_id", all = FALSE)
      colnames(data_sig_info)<-c("GO","pvalues_adjBY", "Term")
      write.table(data_sig_info,file=paste0("topGO_enrichment_", anals, "_", countp, "_sign_descript.txt"), col.names=FALSE,row.names=FALSE,quote=FALSE, sep = "\t")
      # data_sig0.01 <- data[data$pvalues < 0.01]
      # countp0.01 <- nrow(data_sig0.01)
      # data_sig0.01 <- merge(data_sig0.01, GO2term, by.x = "V1", by.y = "go_id", all = FALSE)
      # colnames(data_sig0.01)<-c("GO","pvalues_adjBY", "Term")
      # write.table(data_sig0.01,file=paste0("topGO_enrichment_", anals, "_pvaluesig0.01_", countp0.01, "sign.txt"), col.names=FALSE,row.names=FALSE,quote=FALSE)
      
      #### OBTAINING NICE TREEMAP: code wrote by Gemma I. Martínez-Redondo on February, 2024 
      
      #Obtain significative GOs
      userData <- data_sig
      colnames(userData)<-c("GO","pvalue")
      userData$description<-get_names(userData$GO)[2] %>% unlist()
      head(userData)
      go_ids<-userData$GO
      
      #Clusters generation
      set.seed(1234)
      mat <- GO_similarity(go_ids,ont=ont,measure="Wang") #Compute semantic similarity from scratch
      df<-cluster_terms(mat, method = "kmeans")
      goClusters<-data.frame(colnames(mat),df) #Format GO Cluster
      colnames(goClusters)<-c("GO","Cluster")
      clusters_representatives<-data.frame()
      for (cluster in unique(goClusters$Cluster)){
        cluster_gos<-goClusters[goClusters$Cluster==cluster,1] #Get GO terms belonging to each cluster
        if (length(cluster_gos)==1){
          clusters_representatives<-rbind(clusters_representatives,c(cluster,cluster_gos,get_names(cluster_gos)$go_name))
          next
        }
        submat<-mat[cluster_gos,cluster_gos]
        closeness_values<-graph_from_adjacency_matrix(as.matrix(submat),weighted=TRUE) %>% closeness()
        if(sum(is.na(closeness_values))!=0){
          i<-1
          for (go in names(closeness_values)){
            clusters_representatives<-rbind(clusters_representatives,c(paste0(cluster,"_",i),go,get_names(go)$go_name))
            i<-i+1
          }
        }
        else{
          representative<-closeness_values[which.min(closeness_values)] %>% names() #Obtain "furthest" node, which means that it's the more general one
          clusters_representatives<-rbind(clusters_representatives,c(cluster,representative,get_names(representative)$go_name))
        }
      }
      colnames(clusters_representatives)<-c("Cluster","GO","Description")
      clusters_representatives$Cluster<-as.integer(clusters_representatives$Cluster)
      goClusters<-inner_join(goClusters,clusters_representatives,by="Cluster")
      colnames(goClusters)<-c("GO","Cluster","GO_representative","representative")
      
      #Save information on clusters 
      write.table(goClusters,file=paste0("topGO_enrichment_", anals, "_", countp, "_GO_clusters.txt"), row.names=FALSE,quote=FALSE,col.names=FALSE,sep = "\t")
      
      #Joining cluster's information
      df<-inner_join(userData,goClusters,by="GO")
      df<-df %>% mutate(value=case_when(pvalue<1e-300 ~ 300, pvalue>=1e-300 ~ abs(log10(pvalue)))) %>% select(c(description,representative,value))
      df$value<-as.numeric(as.character(df$value))
      df<-df[,c("description","value","representative")]
      df$representative<-as.factor(df$representative)
      
      #Delete NA from df
      newdf <- df[complete.cases(df), ]
      , you can select the palette that you want
      #Color palette
      color_palette <- scico(length(unique(df$representative)), palette = 'navia')  #green palette, you can select the palette that you want
      plot_title<-paste0("TreeMap_", anals)
      
      # by default, outputs to a PDF file
      pdf( file=paste0(anals, "_TreeMap.pdf"), width=16, height=9 ) # width and height are in inches
      # check the tmPlot command documentation for all possible parameters - there are a lot more
      treemap(
        newdf,
        index = c("representative","description"),
        vSize = "value",
        type = "categorical",
        vColor = "representative",
        title = plot_title,
        inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
        lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
        bg.labels = "#F0EAF9",   # define background color of group labels
        # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
        position.legend = "none",
        palette = color_palette,
        border.col = "white",
        border.lwds = c(4,2)
      )
      dev.off()
      
    }
  }
}