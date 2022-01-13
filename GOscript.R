library(topGO)
set1 <- read.table("~/geneids.txt", h=F) ### Import textfile of gene IDs identified during DE analyses (see DE script)
colnames(set1) <- c("ID","Gene.name")


main.set <- set1$ID[ set1$Gene.name %in% row.names(dds)]


canada <- readMappings(file = "./geneidgo.db") # Import mapped GOterm (from the Salvelinus genome assembly ASM291031v2)




#### targets of differentially expressed miRNAs

mi.de.list <- readRDS("~/de.list.miRNA.RDS") # ID of differentially expressed miRNAs (from the miRNA-de.R script)

mi.de.names<- c("sh-pl-sb-150","sh-plsb-plpl-150", "sh-sbpl-plpl-150" ,"sh-plsb-sbsb-150",
             "sh-sbpl-sbsb-150", "sh-pl-sb-200","sh-plsb-plpl-200", "sh-sbpl-plpl-200", 
             "sh-plsb-sbsb-200", "sh-sbpl-sbsb-200")

GOlist.de.mi <- list()

for(i in 1:length(mi.de.list)){
  sel.set <-df.target$GENEID[df.target$miRNA %in% rownames(na.omit(mi.de.list[[i]])
                                                [na.omit(mi.de.list[[i]]$padj) < 0.1,])]
  sel.set <- sel.set[!duplicated(sel.set)]
  sel.set <- set1$ID[set1$Gene.name %in% sel.set]
  geneList <- factor(as.integer(main.set %in% sel.set))
  names(geneList) <- main.set
  if(length(levels(as.factor(geneList))) == 2) {
    tGOdata <- new("topGOdata",
                   description = mi.de.names[i],
                   ontology = "BP",
                   allGenes = geneList, 
                   nodeSize = 10,
                   annot = annFUN.gene2GO, gene2GO = canada)
    wtest <- runTest(tGOdata, algorithm = "weight01", statistic = "fisher")
    
    gt.de1 <- GenTable(tGOdata, weight01Fisher = wtest, 
                       orderBy = "weight01Fisher", 
                       ranksOf = "weight01Fisher",topNodes = 10)
    GOlist.de.mi[[i]] <- gt.de1
  } else {
    GOlist.de.mi[[i]] <- matrix("NA", ncol = 7, nrow = 1) }
}


GOlist.de.mi[[2]] <- as.data.frame(matrix(NA, nrow = 10,ncol = 7))
colnames(GOlist.de.mi[[2]]) <- colnames(GOlist.de.mi[[1]])
GOlist.de.mi[[10]] <- as.data.frame(matrix(NA, nrow = 10,ncol = 7))
colnames(GOlist.de.mi[[10]]) <- colnames(GOlist.de.mi[[1]])

for(i in 1:length(mi.de.names)){
  GOlist.de.mi[[i]]$contrast <- rep(mi.de.names[i], nrow(GOlist.de.mi[[i]])) 
}

dfGOde.mi <- do.call(rbind,GOlist.de.mi)
write.table(dfGOde.mi, "~/go-DE-miRNA.txt")



##### GO of miRNA LCV clusters

dfcl.mi <- readRDS( "~/Desktop/dfcl.mi.rds") # list of micro RNA clusters (from the miRNA-de.R script)

GOlist.lcv.mi <- list()

for(i in 1:length(levels(as.factor(dfcl.mi$miRNACluster)))){
  sel.set <-dfcl.mi$gene[dfcl.mi$miRNACluster == i ]
  sel.set <-df.target$GENEID[df.target$miRNA %in% sel.set]
  sel.set <- sel.set[!duplicated(sel.set)]
  sel.set <- set1$ID[set1$Gene.name %in% sel.set]
  geneList <- factor(as.integer(main.set %in% sel.set))
  names(geneList) <- main.set
  if(length(levels(as.factor(geneList))) == 2){ tGOdata <- new("topGOdata",
                                                               description = mi.de.names[i],
                                                               ontology = "BP",
                                                               allGenes = geneList, 
                                                               nodeSize = 10,
                                                               annot = annFUN.gene2GO, gene2GO = canada)
  wtest <- runTest(tGOdata, algorithm = "weight01", statistic = "fisher")
  
  gt.de1 <- GenTable(tGOdata, weight01Fisher = wtest, 
                     orderBy = "weight01Fisher", 
                     ranksOf = "weight01Fisher",topNodes = 10)
  GOlist.lcv.mi[[i]] <- gt.de1
  } else {
    GOlist.lcv.mi[[i]] <- matrix("NA", ncol = 7, nrow = 1)
  } 
}

for(i in 1:length(mi.de.names)){
  GOlist.lcv.mi[[i]]$contrast <- rep(mi.de.names[i], nrow(GOlist.lcv.mi[[i]])) 
}

names(GOlist.lcv.mi) <- c(1:10)

for(i in 1:length(names(GOlist.lcv.mi))){
  GOlist.lcv.mi[[i]]$contrast <- rep(names(GOlist.lcv.mi)[i], nrow(GOlist.lcv.mi[[i]])) 
}

dfGO.lcv.mi <- do.call(rbind,GOlist.lcv.mi)

write.table( dfGO.lcv.mi, "~/go-lcv-miRNA.txt")




