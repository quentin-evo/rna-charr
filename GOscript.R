library(topGO)
set1 <- read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/DE-analyses/Go-analysis/geneids.txt", h=F)
colnames(set1) <- c("ID","Gene.name")


main.set <- set1$ID[ set1$Gene.name %in% row.names(dds)]


canada <- readMappings(file = "./Desktop/Databases_GO/geneidgo.db") # Import mapped GOterm


##### GO of differentially expressed genes at 200ts

GOlist.de <- list()
de.genes <- c(sh.pl.sb.200, sh.plsb.plpl.200, sh.sbpl.plpl.200, sh.plsb.sbsb.200,
              sh.sbpl.sbsb.200, sh.plsb.sbpl.200)
de.nm <- c("SBxSB/PLxPL-200ts", "PLxSB/PLxPL-200ts", "SBxPL/PLxPL-200ts", "PLxSB-SBxSB-200ts",
              "SBxPL/SBxSV-200ts", "PLxSB/SBxPL-200ts")


# Produce GO list (the object "contr" is from the DE.R manuscript)
for(i in 1:length(names(contr))){
  sel.set <-set1$ID[set1$Gene.name %in% rownames(na.omit(contr[[i]])
                                                 [na.omit(contr[[i]]$padj) < 0.1,])]
  
  geneList <- factor(as.integer(main.set %in% sel.set))
  names(geneList) <- main.set
  
  tGOdata <- new("topGOdata",
                 description = names(contr)[i],
                 ontology = "BP",
                 allGenes = geneList, 
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = canada)
  wtest <- runTest(tGOdata, algorithm = "weight01", statistic = "fisher")
  
  gt.de1 <- GenTable(tGOdata, weight01Fisher = wtest, 
                     orderBy = "weight01Fisher", 
                     ranksOf = "weight01Fisher",topNodes = 10)
  GOlist.de[[i]] <- gt.de1
}

for(i in 1:length(names(contr))){
  GOlist.de[[i]]$contrast <- rep(names(contr)[i], nrow(GOlist.de[[i]])) 
}

dfGO <- do.call(rbind,GOlist.de)

write.table(dfGO, "~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/DE-analyses/Go-analysis/go-DE.txt")


showSigOfNodes(tGOdata, score(ftest), firstSigNodes = 5, useInfo = 'all')
showSigOfNodes(tGOdata, score(ftest2), firstSigNodes = 5, useInfo = 'all')



#### targets of differentially expressed miRNAs

mi.de.list <- readRDS("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/miRNA-target/de.list.miRNA.RDS")

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
write.table(dfGOde.mi, "~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/DE-analyses/Go-analysis/go-DE-miRNA.txt")



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

write.table( dfGO.lcv.mi, "~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/DE-analyses/Go-analysis/go-lcv-miRNA.txt")

#### GO of mRNA LCV clusters

GOlist.lcv.m <- list()

for(i in 1:length(levels(as.factor(dfcl$mRNACluster)))){
  sel.genes <-dfcl$gene[dfcl$mRNACluster == i ]
  sel.genes <- sel.genes[!duplicated(sel.genes)]
  sel.set <-set1$ID[set1$Gene.name %in% sel.genes]
  geneList <- factor(as.integer(main.set %in% sel.set))
  names(geneList) <- main.set
  if(length(levels(as.factor(geneList))) == 2) {
    tGOdata <- new("topGOdata",
                   description = paste("Cluster",i),
                   ontology = "BP",
                   allGenes = geneList, 
                   nodeSize = 10,
                   annot = annFUN.gene2GO, gene2GO = canada)
    wtest <- runTest(tGOdata, algorithm = "weight01", statistic = "fisher")
    
    gt.de1 <- GenTable(tGOdata, weight01Fisher = wtest, 
                       orderBy = "weight01Fisher", 
                       ranksOf = "weight01Fisher",topNodes = 10)
    GOlist.lcv.m[[i]] <- gt.de1
  } else {
    GOlist.lcv.m[[i]] <- matrix("NA", ncol = 7, nrow = 1) }
}
names(GOlist.lcv.m) <- c(1:10)

for(i in 1:length(names(GOlist.lcv.m))){
  GOlist.lcv.m[[i]]$contrast <- rep(names(GOlist.lcv.m)[i], nrow(GOlist.lcv.m[[i]])) 
}

dfGO.lcv.m <- do.call(rbind,GOlist.lcv.m)
write.table(dfGO.lcv.m, "~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/DE-analyses/Go-analysis/go-lvc.mRNA.txt")


