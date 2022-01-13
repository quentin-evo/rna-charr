library(tidyr)
library("DESeq2")
library("tximport")
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(edgeR)
library(ComplexHeatmap)
library(circlize)
library(MCMCglmm)
library(gridExtra)

load("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/DE-analyses/DE-mRNA.RData")

# Import read count
mdf<-read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mimir-saved/mrna_featureCounts.txt",h=T,skip = 1)
rownames(mdf) <- mdf$Geneid
mdf<-mdf[,-c(1,2,3,4,5,6)] # Remove Chr, strand info, lengths
colnames(mdf) <- gsub('X.home.qjb1.rnaseq.STAR.AlignedSortedBAM.', '', colnames(mdf))
colnames(mdf) <- gsub('.STAR.alignedAligned.sortedByCoord.out.bam', '', colnames(mdf))
colnames(mdf) <- gsub('[.]', '_', colnames(mdf))

# Import metadata
df.info<-read.csv("/Users/Muscardin/Documents/Analyses/Post-zygotic mechanisms/RNAseq/hybridsmiRNAseq/hybridsmiRNAseq/dataInfo/hybridsCharrDataFrameDagny.csv",h=T, sep = ";")
df.info$Batch<-as.factor(df.info$Batch)
df.info$Timepoint<-as.factor(df.info$Timepoint)
row.names(df.info)<-df.info$Sample
df.info<-df.info[,2:6]
df.info$group<-paste(df.info$Cross,df.info$Timepoint,sep = "_")
df.info<-df.info[colnames(mdf),]

stopifnot(colnames(mdf) == rownames(df.info))

# Generate DeSeq data set
dds <- DESeqDataSetFromMatrix(countData = mdf,
                              colData = df.info,
                              design = ~Batch + group)

# Discard genes with less than 10 read counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

##### Exploratory analyses 
vsd<-vst(dds, blind=FALSE)  # Transformation by Variance Stabilizing Transformation

# PCA 
pcaData <- plotPCA(vsd, intgroup=c("group", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$group.1<-factor(pcaData$group.1, levels = c("PLxPL_150","SBxSB_150","PLxSB_150","SBxPL_150","PLxPL_200","SBxSB_200","PLxSB_200","SBxPL_200"))
ggplot(pcaData, aes(PC1, PC2, color=group.1, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("#003300", "#336600", "99FF99","#99FFCC",
                                "#003366", "#0066CC", "#99CCFF", "#CCE5FF")) +
  labs(color= "Cross & Age") +
  theme_bw()

pcaData[pcaData$PC2 < -15,] # "Extreme" samples along PC2
write.csv("~/Desktop/PCA-samples.txt",x = pcaData[pcaData$PC2 < -15,])


# Heatmap
select <- order(rowMeans(counts(dds,normalized=F)),
                decreasing=TRUE)[1:100]
# Specify colors
ann_colors = list(
  group =c(PLxPL_150 = "#003300", SBxSB_150 = "#336600", PLxSB_150 = "#99FF99",SBxPL_150 = "#99FFCC",
           PLxPL_200 = "#003366", SBxSB_200 ="#0066CC", PLxSB_200 ="#99CCFF", SBxPL_200 = "#CCE5FF"),
  Batch = c("1" = "#202020", "2" = "#808080", "3" = "#E0E0E0")
)
hdf <- as.data.frame(colData(dds)[,c("group","Batch")])
pheatmap(assay(vsd), cluster_rows=T, show_rownames=FALSE,
         cluster_cols=T, annotation_col=hdf ,
         annotation_colors = ann_colors
) 
# 149_200_C mixed up with the 150ts samples. In the PCA, it has coordinates: 2.564065 ; -19.8883017


# DE analysis
ds_de<-DESeq(dds)
resultsNames(ds_de)

# 150ts
r.pl.sb.150 <- results(ds_de,contrast =  c("group","SBxSB_150","PLxPL_150"))
r.plsb.plpl.150 <- results(ds_de,contrast =  c("group","PLxSB_150","PLxPL_150"))
r.sbpl.plpl.150 <- results(ds_de,contrast =  c("group","SBxPL_150","PLxPL_150"))
r.plsb.sbsb.150 <- results(ds_de,contrast =  c("group","PLxSB_150","SBxSB_150"))
r.sbpl.sbsb.150 <- results(ds_de,contrast =  c("group","SBxPL_150","SBxSB_150"))
r.sbpl.sbsb.150 <- results(ds_de,contrast =  c("group","SBxPL_150","SBxSB_150"))
r.plsb.sbpl.150 <- results(ds_de,contrast =  c("group","PLxSB_150","SBxPL_150"))

# 200ts
r.pl.sb.200 <- results(ds_de,contrast =  c("group","SBxSB_200","PLxPL_200"))
r.plsb.plpl.200 <- results(ds_de,contrast =  c("group","PLxSB_200","PLxPL_200"))
r.sbpl.plpl.200 <- results(ds_de,contrast =  c("group","SBxPL_200","PLxPL_200"))
r.plsb.sbsb.200 <- results(ds_de,contrast =  c("group","PLxSB_200","SBxSB_200"))
r.sbpl.sbsb.200 <- results(ds_de,contrast =  c("group","SBxPL_200","SBxSB_200"))
r.plsb.sbpl.200 <- results(ds_de,contrast =  c("group","PLxSB_200","SBxPL_200"))

# Apply shrinkage
dds$group<-relevel(dds$group,"SBxSB_200")
ds.de.2<-nbinomWaldTest(ds_de)

# 150ts
contr<-list()

contr[["sh.pl.sb.150"]] <- sh.pl.sb.150.2 <- sh.pl.sb.150.2 <- sh.pl.sb.150 <- lfcShrink(ds_de, type="ashr", res = r.pl.sb.150, parallel = T,BPPARAM=MulticoreParam(4))
contr[["sh.plsb.plpl.150"]] <- sh.plsb.plpl.150 <- lfcShrink(ds_de, type="ashr", res = r.plsb.plpl.150, parallel = T,BPPARAM=MulticoreParam(4))
contr[["sh.sbpl.plpl.150"]]  <- sh.sbpl.plpl.150 <- lfcShrink(ds_de, type="ashr", res = r.sbpl.plpl.150, parallel = T,BPPARAM=MulticoreParam(4))
contr[["sh.plsb.sbsb.150"]]  <- sh.plsb.sbsb.150<- lfcShrink(ds_de, type="ashr", res = r.plsb.sbsb.150, parallel = T,BPPARAM=MulticoreParam(4))
contr[["sh.sbpl.sbsb.150"]]  <- sh.sbpl.sbsb.150<- lfcShrink(ds_de, type="ashr", res = r.sbpl.sbsb.150, parallel = T,BPPARAM=MulticoreParam(4))
contr[["sh.plsb.sbpl.150"]]  <- sh.plsb.sbpl.150<- lfcShrink(ds_de, type="ashr", res = r.plsb.sbpl.150, parallel = T,BPPARAM=MulticoreParam(4))

# 200ts
contr[["sh.pl.sb.200"]]  <- sh.pl.sb.200 <- lfcShrink(ds_de, type="ashr", res = r.pl.sb.200, parallel = T,BPPARAM=MulticoreParam(4))
contr[["sh.plsb.plpl.200"]]  <- sh.plsb.plpl.200<- lfcShrink(ds_de, type="ashr", res = r.plsb.plpl.200, parallel = T,BPPARAM=MulticoreParam(4))
contr[["sh.sbpl.plpl.200"]]  <-sh.sbpl.plpl.200  <- lfcShrink(ds_de, type="ashr", res = r.sbpl.plpl.200, parallel = T,BPPARAM=MulticoreParam(4))
contr[["sh.plsb.sbsb.200"]]  <- sh.plsb.sbsb.200<- lfcShrink(ds_de, type="ashr", res = r.plsb.sbsb.200, parallel = T,BPPARAM=MulticoreParam(4))
contr[["sh.sbpl.sbsb.200"]]  <- sh.sbpl.sbsb.200<- lfcShrink(ds_de, type="ashr", res = r.sbpl.sbsb.200, parallel = T,BPPARAM=MulticoreParam(4))
contr[["sh.plsb.sbpl.200"]]  <- sh.plsb.sbpl.200<- lfcShrink(ds_de, type="ashr", res = r.plsb.sbpl.200, parallel = T,BPPARAM=MulticoreParam(4))

##### MA plot ################

# PLxPL vs SBxSB 
par(mfrow = c(1,2))
DESeq2::plotMA(sh.pl.sb.150,  ylim=c(-5.5,5.5))
title(main="a.      SBxSB/PLxPL (150ts)"
      , adj = 0
)
DESeq2::plotMA(sh.pl.sb.200,ylim=c(-5.5,5.5))
title(main="b.      SBxSB/PLxPL (200ts)"
      , adj = 0,
)

# hybrids 150ts
par(mfrow = c(2,2))
DESeq2::plotMA(sh.plsb.plpl.150, ylim=c(-5,5))
title(main="a.  PLxSB/PLxPL (150ts)", adj = 0)
DESeq2::plotMA(sh.plsb.sbsb.150, ylim=c(-5,5))
title(main="c.  PLxSB/SBxSB (150ts)", adj = 0)
DESeq2::plotMA(sh.sbpl.plpl.150, ylim=c(-5,5))
title(main="b.  SBxPL/PLxPL (150ts)", adj = 0)
DESeq2::plotMA(sh.sbpl.sbsb.150, ylim=c(-5,5))
title(main="d.  SBxPL/SBxSB (150ts)", adj = 0)

# hybrids 200ts
par(mfrow = c(2,2))
DESeq2::plotMA(sh.plsb.plpl.200, ylim=c(-5,5))
title(main="a.  PLxSB/PLxPL (200ts)", adj = 0)
DESeq2::plotMA(sh.plsb.sbsb.200, ylim=c(-5,5))
title(main="c.  PLxSB/SBxSB (200ts)", adj = 0)
DESeq2::plotMA(sh.sbpl.plpl.200, ylim=c(-5,5))
title(main="b.  SBxPL/PLxPL (200ts)", adj = 0)
DESeq2::plotMA(sh.sbpl.sbsb.200, ylim=c(-5,5))
title(main="d.  SBxPL/SBxSB (200ts)", adj = 0)

###### Maternal inheritance test ##########

fmatern <- function(pop1, pop2,name1,name2,ncomp) {
  cont<- matrix(ncol = 2,nrow = 2)
  cont[1,1] <- sum(pop1$padj < 0.1,na.rm = T)
  cont[2,1] <- sum(pop1$baseMean > 0,na.rm = T) - cont[1,1]
  cont[1,2] <- sum(pop2$padj < 0.1,na.rm = T) 
  cont[2,2] <- sum(pop2$baseMean > 0,na.rm = T) - cont[1,2]
  chisq <- chisq.test(cont)
  return(data.frame(
    contrast1 = paste(name1 , " = ",cont[1,1], " (", round((cont[1,1]/cont[2,1])*100,2),"%)",sep=""),
    contrast2 = paste(name2 , " = ",cont[1,2], " (", round((cont[1,2]/cont[2,2])*100,2),"%)",sep=""),
    Xsquared = chisq$statistic,
    p = round(chisq$p.value,2),
    Bonferroni = ifelse(test = chisq$p.value < 0.05/ncomp,yes = "*",no = " ")))
}

chitable <- rbind(
  fmatern(pop1 = sh.plsb.plpl.150, pop2 = sh.sbpl.plpl.150, 
        name1 = "PLxSB/PLxPL-150ts", name2 = "SBxPL/PLxPL-150ts",ncomp=4),
  fmatern(pop1 = sh.plsb.sbsb.150, pop2 = sh.sbpl.sbsb.150, 
          name1 = "PLxSB/SBxSB-150ts", name2 = "SBxPL/SBxSB-150ts",ncomp=4),
  fmatern(pop1 = sh.plsb.plpl.200, pop2 = sh.sbpl.plpl.200, 
          name1 = "PLxSB/PLxPL-200ts", name2 = "SBxPL/PLxPL-200ts",ncomp=4),
  fmatern(pop1 = sh.plsb.sbsb.200, pop2 = sh.sbpl.sbsb.200, 
          name1 = "PLxSB/SBxSB-200ts", name2 = "SBxPL/SBxSB-200ts",ncomp=4))
write.csv("~/Desktop/maternal-test.txt", x = chitable )

contPLPL<- matrix(ncol = 2,nrow = 2)
contPLPL[1,1] <- sum(sh.plsb.plpl.150$padj < 0.1,na.rm = T)
contPLPL[2,1] <- sum(sh.plsb.plpl.150$baseMean > 0,na.rm = T) - contPLPL[1,1]
contPLPL[1,2] <- sum(sh.sbpl.plpl.150$padj < 0.1,na.rm = T) 
contPLPL[2,2] <- sum(sh.sbpl.plpl.150$baseMean > 0,na.rm = T) - contPLPL[1,2]
chiPL200 <- chisq.test(contPLPL)

contPLPL[1,1] <- sum(sh.pl.sb.150$padj < 0.1,na.rm = T)/sum(sh.pl.sb.150$baseMean > 0,na.rm = T)

#### Inheritance patterns 

inh <- function(df) {
  dfcontr <- as.data.frame(df)
  dfcontr$inheritance150 <- rep(NA, nrow(dfcontr))
  dfcontr$inheritance200 <- rep(NA, nrow(dfcontr))
  dfcontr$inheritance150 <- as.factor(dfcontr$inheritance150)
  levels(dfcontr$inheritance150) = c("Maternal","PLdominant", "SBdominant", "Additive",
                             "Overdominant","Underdominant","Unclassified","NoDE")
  dfcontr$inheritance200 <- as.factor(dfcontr$inheritance200)
  levels(dfcontr$inheritance200) = c("Maternal","PLdominant", "SBdominant", "Additive",
                             "Overdominant","Underdominant","Unclassified","NoDE")
  ## 150 ts
  # Maternal

  dfcontr$inheritance150[dfcontr$sh.pl.sb.150.padj < 0.1
                         & dfcontr$sh.plsb.plpl.150.padj > 0.1
                         & dfcontr$sh.sbpl.plpl.150.padj < 0.1
                         & dfcontr$sh.plsb.sbsb.150.padj < 0.1
                         & dfcontr$sh.sbpl.sbsb.150.padj > 0.1
                         & dfcontr$sh.pl.sb.150.log2FoldChange < 0
                         & dfcontr$sh.sbpl.plpl.150.log2FoldChange < 0
                         & dfcontr$sh.plsb.sbsb.150.log2FoldChange > 0
  ] <-"Maternal"
  dfcontr$inheritance150[dfcontr$sh.pl.sb.150.padj < 0.1
                         & dfcontr$sh.plsb.plpl.150.padj > 0.1
                         & dfcontr$sh.sbpl.plpl.150.padj < 0.1
                         & dfcontr$sh.plsb.sbsb.150.padj < 0.1
                         & dfcontr$sh.sbpl.sbsb.150.padj > 0.1
                         & dfcontr$sh.pl.sb.150.log2FoldChange > 0
                         & dfcontr$sh.sbpl.plpl.150.log2FoldChange > 0
                         & dfcontr$sh.plsb.sbsb.150.log2FoldChange < 0
  ] <- "Maternal"

  # PL-dominant
  dfcontr$inheritance150[dfcontr$sh.pl.sb.150.padj < 0.1
                         & dfcontr$sh.plsb.plpl.150.padj > 0.1
                         & dfcontr$sh.sbpl.plpl.150.padj > 0.1
                         & dfcontr$sh.plsb.sbsb.150.padj < 0.1
                         & dfcontr$sh.sbpl.sbsb.150.padj < 0.1
                         & dfcontr$sh.pl.sb.150.log2FoldChange < 0
                         & dfcontr$sh.sbpl.sbsb.150.log2FoldChange > 0
                         & dfcontr$sh.plsb.sbsb.150.log2FoldChange > 0
  ] <-"PLdominant"

  dfcontr$inheritance150[dfcontr$sh.pl.sb.150.padj < 0.1
                         & dfcontr$sh.plsb.plpl.150.padj > 0.1
                         & dfcontr$sh.sbpl.plpl.150.padj > 0.1
                         & dfcontr$sh.plsb.sbsb.150.padj < 0.1
                         & dfcontr$sh.sbpl.sbsb.150.padj < 0.1
                         & dfcontr$sh.pl.sb.150.log2FoldChange > 0
                         & dfcontr$sh.sbpl.sbsb.150.log2FoldChange < 0
                         & dfcontr$sh.plsb.sbsb.150.log2FoldChange < 0
  ] <-"PLdominant"

  dfcontr$inheritance150[dfcontr$sh.pl.sb.150.padj < 0.1
                         & dfcontr$sh.plsb.sbsb.150.padj > 0.1
                         & dfcontr$sh.sbpl.sbsb.150.padj > 0.1
                         & dfcontr$sh.plsb.plpl.150.padj < 0.1
                         & dfcontr$sh.sbpl.plpl.150.padj < 0.1
                         & dfcontr$sh.pl.sb.150.log2FoldChange < 0
                         & dfcontr$sh.sbpl.plpl.150.log2FoldChange > 0
                         & dfcontr$sh.plsb.plpl.150.log2FoldChange > 0
  ] <-"SBdominant"

  dfcontr$inheritance150[dfcontr$sh.pl.sb.150.padj < 0.1
                         & dfcontr$sh.plsb.sbsb.150.padj > 0.1
                         & dfcontr$sh.sbpl.sbsb.150.padj > 0.1
                         & dfcontr$sh.plsb.plpl.150.padj < 0.1
                         & dfcontr$sh.sbpl.plpl.150.padj < 0.1
                         & dfcontr$sh.pl.sb.150.log2FoldChange > 0
                         & dfcontr$sh.sbpl.plpl.150.log2FoldChange < 0
                         & dfcontr$sh.plsb.plpl.150.log2FoldChange < 0
  ] <-"SBdominant"

  # Additive
  dfcontr$inheritance150[dfcontr$sh.pl.sb.150.padj < 0.1
                         & dfcontr$sh.plsb.sbsb.150.padj < 0.1
                         & dfcontr$sh.sbpl.sbsb.150.padj < 0.1
                         & dfcontr$sh.plsb.plpl.150.padj < 0.1
                         & dfcontr$sh.sbpl.plpl.150.padj < 0.1
                         & dfcontr$sh.pl.sb.150.log2FoldChange < 0
                         & dfcontr$sh.sbpl.plpl.150.log2FoldChange < 0
                         & dfcontr$sh.plsb.plpl.150.log2FoldChange < 0
                         & dfcontr$sh.sbpl.sbsb.150.log2FoldChange > 0
                         & dfcontr$sh.plsb.sbsb.150.log2FoldChange > 0
  ] <-"Additive"
  dfcontr$inheritance150[dfcontr$sh.pl.sb.150.padj < 0.1
                         & dfcontr$sh.plsb.sbsb.150.padj < 0.1
                         & dfcontr$sh.sbpl.sbsb.150.padj < 0.1
                         & dfcontr$sh.plsb.plpl.150.padj < 0.1
                         & dfcontr$sh.sbpl.plpl.150.padj < 0.1
                         & dfcontr$sh.pl.sb.150.log2FoldChange > 0
                         & dfcontr$sh.sbpl.plpl.150.log2FoldChange > 0
                         & dfcontr$sh.plsb.plpl.150.log2FoldChange > 0
                         & dfcontr$sh.sbpl.sbsb.150.log2FoldChange < 0
                         & dfcontr$sh.plsb.sbsb.150.log2FoldChange < 0
  ] <-"Additive"

  #Overdominant
  dfcontr$inheritance150[ dfcontr$sh.plsb.sbsb.150.padj < 0.1
                          & dfcontr$sh.sbpl.sbsb.150.padj < 0.1
                          & dfcontr$sh.plsb.plpl.150.padj < 0.1
                          & dfcontr$sh.sbpl.plpl.150.padj < 0.1
                          & dfcontr$sh.sbpl.plpl.150.log2FoldChange > 0
                          & dfcontr$sh.plsb.plpl.150.log2FoldChange > 0
                          & dfcontr$sh.sbpl.sbsb.150.log2FoldChange > 0
                          & dfcontr$sh.plsb.sbsb.150.log2FoldChange > 0
  ] <-"Overdominant"
  # Underdominant
  dfcontr$inheritance150[ dfcontr$sh.plsb.sbsb.150.padj < 0.1
                          & dfcontr$sh.sbpl.sbsb.150.padj < 0.1
                          & dfcontr$sh.plsb.plpl.150.padj < 0.1
                          & dfcontr$sh.sbpl.plpl.150.padj < 0.1
                          & dfcontr$sh.sbpl.plpl.150.log2FoldChange < 0
                          & dfcontr$sh.plsb.plpl.150.log2FoldChange < 0
                          & dfcontr$sh.sbpl.sbsb.150.log2FoldChange < 0
                          & dfcontr$sh.plsb.sbsb.150.log2FoldChange < 0
  ] <-"Underdominant"

  # No DE
  dfcontr$inheritance150[ dfcontr$sh.plsb.sbsb.150.padj > 0.1
                          & dfcontr$sh.plsb.sbsb.150.padj > 0.1
                          & dfcontr$sh.sbpl.sbsb.150.padj > 0.1
                          & dfcontr$sh.plsb.plpl.150.padj > 0.1
                          & dfcontr$sh.sbpl.plpl.150.padj > 0.1
  ] <-"NoDE"

############### 200ts

#dfcontr$inheritance200 <- rep(NA,nrow(dfcontr))

#is.na(dfcontr[ ,c("sh.pl.sb.200.padj","sh.plsb.sbsb.200.padj", "sh.sbpl.sbsb.200.padj","sh.plsb.plpl.200.padj","sh.sbpl.plpl.200.padj")]) <- 1

# Maternal
dfcontr$inheritance200[dfcontr$sh.pl.sb.200.padj < 0.1
                       & dfcontr$sh.plsb.plpl.200.padj > 0.1
                       & dfcontr$sh.sbpl.plpl.200.padj < 0.1
                       & dfcontr$sh.plsb.sbsb.200.padj < 0.1
                       & dfcontr$sh.sbpl.sbsb.200.padj > 0.1
                       & dfcontr$sh.pl.sb.200.log2FoldChange < 0
                       & dfcontr$sh.sbpl.plpl.200.log2FoldChange < 0
                       & dfcontr$sh.plsb.sbsb.200.log2FoldChange > 0
] <-"Maternal"
dfcontr$inheritance200[dfcontr$sh.pl.sb.200.padj < 0.1
                       & dfcontr$sh.plsb.plpl.200.padj > 0.1
                       & dfcontr$sh.sbpl.plpl.200.padj < 0.1
                       & dfcontr$sh.plsb.sbsb.200.padj < 0.1
                       & dfcontr$sh.sbpl.sbsb.200.padj > 0.1
                       & dfcontr$sh.pl.sb.200.log2FoldChange > 0
                       & dfcontr$sh.sbpl.plpl.200.log2FoldChange > 0
                       & dfcontr$sh.plsb.sbsb.200.log2FoldChange < 0
] <- "Maternal"

# PL-dominant
dfcontr$inheritance200[is.na(dfcontr$inheritance200) 
                       & dfcontr$sh.pl.sb.200.padj < 0.1
                       & dfcontr$sh.plsb.plpl.200.padj > 0.1
                       & dfcontr$sh.sbpl.plpl.200.padj > 0.1
                       & dfcontr$sh.plsb.sbsb.200.padj < 0.1
                       & dfcontr$sh.sbpl.sbsb.200.padj < 0.1
                       & dfcontr$sh.pl.sb.200.log2FoldChange < 0
                       & dfcontr$sh.sbpl.sbsb.200.log2FoldChange > 0
                       & dfcontr$sh.plsb.sbsb.200.log2FoldChange > 0
] <-"PLdominant"

dfcontr$inheritance200[is.na(dfcontr$inheritance200) 
                       & dfcontr$sh.pl.sb.200.padj < 0.1
                       & dfcontr$sh.plsb.plpl.200.padj > 0.1
                       & dfcontr$sh.sbpl.plpl.200.padj > 0.1
                       & dfcontr$sh.plsb.sbsb.200.padj < 0.1
                       & dfcontr$sh.sbpl.sbsb.200.padj < 0.1
                       & dfcontr$sh.pl.sb.200.log2FoldChange > 0
                       & dfcontr$sh.sbpl.sbsb.200.log2FoldChange < 0
                       & dfcontr$sh.plsb.sbsb.200.log2FoldChange < 0
] <-"PLdominant"

dfcontr$inheritance200[dfcontr$sh.pl.sb.200.padj < 0.1
                       & dfcontr$sh.plsb.sbsb.200.padj > 0.1
                       & dfcontr$sh.sbpl.sbsb.200.padj > 0.1
                       & dfcontr$sh.plsb.plpl.200.padj < 0.1
                       & dfcontr$sh.sbpl.plpl.200.padj < 0.1
                       & dfcontr$sh.pl.sb.200.log2FoldChange < 0
                       & dfcontr$sh.sbpl.plpl.200.log2FoldChange > 0
                       & dfcontr$sh.plsb.plpl.200.log2FoldChange > 0
] <-"SBdominant"

dfcontr$inheritance200[dfcontr$sh.pl.sb.200.padj < 0.1
                       & dfcontr$sh.plsb.sbsb.200.padj > 0.1
                       & dfcontr$sh.sbpl.sbsb.200.padj > 0.1
                       & dfcontr$sh.plsb.plpl.200.padj < 0.1
                       & dfcontr$sh.sbpl.plpl.200.padj < 0.1
                       & dfcontr$sh.pl.sb.200.log2FoldChange > 0
                       & dfcontr$sh.sbpl.plpl.200.log2FoldChange < 0
                       & dfcontr$sh.plsb.plpl.200.log2FoldChange < 0
] <-"SBdominant"

# Additive
dfcontr$inheritance200[dfcontr$sh.pl.sb.200.padj < 0.1
                       & dfcontr$sh.plsb.sbsb.200.padj < 0.1
                       & dfcontr$sh.sbpl.sbsb.200.padj < 0.1
                       & dfcontr$sh.plsb.plpl.200.padj < 0.1
                       & dfcontr$sh.sbpl.plpl.200.padj < 0.1
                       & dfcontr$sh.pl.sb.200.log2FoldChange < 0
                       & dfcontr$sh.sbpl.plpl.200.log2FoldChange < 0
                       & dfcontr$sh.plsb.plpl.200.log2FoldChange < 0
                       & dfcontr$sh.sbpl.sbsb.200.log2FoldChange > 0
                       & dfcontr$sh.plsb.sbsb.200.log2FoldChange > 0
] <-"Additive"
dfcontr$inheritance200[dfcontr$sh.pl.sb.200.padj < 0.1
                       & dfcontr$sh.plsb.sbsb.200.padj < 0.1
                       & dfcontr$sh.sbpl.sbsb.200.padj < 0.1
                       & dfcontr$sh.plsb.plpl.200.padj < 0.1
                       & dfcontr$sh.sbpl.plpl.200.padj < 0.1
                       & dfcontr$sh.pl.sb.200.log2FoldChange > 0
                       & dfcontr$sh.sbpl.plpl.200.log2FoldChange > 0
                       & dfcontr$sh.plsb.plpl.200.log2FoldChange > 0
                       & dfcontr$sh.sbpl.sbsb.200.log2FoldChange < 0
                       & dfcontr$sh.plsb.sbsb.200.log2FoldChange < 0
] <-"Additive"

#Overdominant
dfcontr$inheritance200[ dfcontr$sh.plsb.sbsb.200.padj < 0.1
                        & dfcontr$sh.sbpl.sbsb.200.padj < 0.1
                        & dfcontr$sh.plsb.plpl.200.padj < 0.1
                        & dfcontr$sh.sbpl.plpl.200.padj < 0.1
                        & dfcontr$sh.sbpl.plpl.200.log2FoldChange > 0
                        & dfcontr$sh.plsb.plpl.200.log2FoldChange > 0
                        & dfcontr$sh.sbpl.sbsb.200.log2FoldChange > 0
                        & dfcontr$sh.plsb.sbsb.200.log2FoldChange > 0
] <-"Overdominant"
# Underdominant
dfcontr$inheritance200[ dfcontr$sh.plsb.sbsb.200.padj < 0.1
                        & dfcontr$sh.sbpl.sbsb.200.padj < 0.1
                        & dfcontr$sh.plsb.plpl.200.padj < 0.1
                        & dfcontr$sh.sbpl.plpl.200.padj < 0.1
                        & dfcontr$sh.sbpl.plpl.200.log2FoldChange < 0
                        & dfcontr$sh.plsb.plpl.200.log2FoldChange < 0
                        & dfcontr$sh.sbpl.sbsb.200.log2FoldChange < 0
                        & dfcontr$sh.plsb.sbsb.200.log2FoldChange < 0
] <-"Underdominant"

# # No DE
# dfcontr$inheritance200[ dfcontr$sh.plsb.sbsb.200.padj > 0.1
#                         & dfcontr$sh.plsb.sbsb.200.padj > 0.1
#                         & dfcontr$sh.sbpl.sbsb.200.padj > 0.1
#                         & dfcontr$sh.plsb.plpl.200.padj > 0.1
#                         & dfcontr$sh.sbpl.plpl.200.padj > 0.1
# ] <-"NoDE"
# dfcontr$inheritance200[ is.na(dfcontr$sh.plsb.sbsb.200.padj)
#                         | is.na(dfcontr$sh.plsb.sbsb.200.padj)
#                         | is.na(dfcontr$sh.sbpl.sbsb.200.padj)
#                         | is.na( dfcontr$sh.plsb.plpl.200.padj)
#                         | is.na(dfcontr$sh.sbpl.plpl.200.padj)
# ] <-"Unclassified"
# dfcontr$inheritance150 <- as.factor(dfcontr$inheritance150)
# levels(dfcontr$inheritance150) = c("Maternal","PLdominant", "SBdominant", "Additive",
# "Overdominant","Underdominant","Unclassified","NoDE")
# dfcontr$inheritance200 <- as.factor(dfcontr$inheritance150)
# levels(dfcontr$inheritance200) = c("Maternal","PLdominant", "SBdominant", "Additive",
#                            "Overdominant","Underdominant","Unclassified","NoDE")
#return(data.frame("dfcontr$inheritance" = levels(dfcontr$inheritance150), "ts150ts" = summary( as.factor(dfcontr$inheritance150)),
#                  "ts200" = summary(as.factor(dfcontr$inheritance150))))
return(dfcontr)
}
  
dfcontr <- inh(contr)
summary(dfcontr$inheritance200)

### Differentially expressed genes among morphs
write("~/Desktop/Go-analysis/reference-list.txt",x = rownames(sh.pl.sb.150))
write("~/Desktop/Go-analysis/GO-gene-list-target-sbpl-150.txt",x = rownames(na.omit(sh.pl.sb.150)[na.omit(sh.pl.sb.150$padj) < 0.1,]))
write("~/Desktop/Go-analysis/GO-gene-list-target-sbpl-200.txt",x = rownames(na.omit(sh.pl.sb.200)[na.omit(sh.pl.sb.200$padj) < 0.1,]))


### Differentially expressed genes among morphs and hybrids
write("~/Desktop/GO-gene-list-sbpl_plpl-150.txt",x = rownames(sh.sbpl.plpl.150[order(sh.sbpl.plpl.150$padj),]))
write("~/Desktop/Go-analysis/GO-gene-list-target-sbpl_plpl-150.txt",x = rownames(na.omit(sh.sbpl.plpl.150)[na.omit(sh.sbpl.plpl.150$padj) < 0.1,]))
write("~/Desktop/Go-analysis/GO-gene-list-target-sbpl_plpl-200.txt",x = rownames(na.omit(sh.sbpl.plpl.200)[na.omit(sh.sbpl.plpl.200$padj) < 0.1,]))
write("~/Desktop/Go-analysis/GO-gene-list-target-plsb_plpl-200.txt",x = rownames(na.omit(sh.plsb.plpl.200)[na.omit(sh.plsb.plpl.200$padj) < 0.1,]))

write("~/Desktop/Go-analysis/GO-gene-list-target-sbpl_sbsb-200.txt",x = rownames(na.omit(sh.sbpl.sbsb.200)[na.omit(sh.sbpl.sbsb.200$padj) < 0.1,]))
write("~/Desktop/Go-analysis/GO-gene-list-target-plsb_sbsb-200.txt",x = rownames(na.omit(sh.plsb.sbsb.200)[na.omit(sh.plsb.sbsb.200$padj) < 0.1,]))

####### Select Maternaly inherited genes
sel <- sh.pl.sb.200[rownames(dfcontr[complete.cases(dfcontr$inheritance200) & dfcontr$inheritance200 == "Maternal",]),] 
DESeq2::plotMA(sel, 
               ylim=c(-5,5), xlim =c(1e-02,1e+05))
idx <- identify(sel$baseMean, 
                sel$log2FoldChange) # Indentify candidate miRNA interractively
rownames(sel)[idx]

rownames(rbind(sh.pl.sb.200[rownames(dfcontr[complete.cases(dfcontr$inheritance200) &
                                               dfcontr$inheritance200 == "Maternal" & dfcontr$sh.pl.sb.200.log2FoldChange > 1,]),],
      sh.pl.sb.200[rownames(dfcontr[complete.cases(dfcontr$inheritance200) &
                                      dfcontr$inheritance200 == "Maternal" & dfcontr$sh.pl.sb.200.log2FoldChange < -1,]),]))

write("~/Desktop/Go-analysis/maternal200.txt",x = rownames(sel))


####### Select PL-dominant genes
selpl <- sh.pl.sb.200[rownames(dfcontr[complete.cases(dfcontr$inheritance200) & dfcontr$inheritance200 == "PLdominant",]),] 
DESeq2::plotMA(selpl, 
               ylim=c(-5,5), xlim =c(1e-02,1e+05))
idx <- identify(selpl$baseMean, 
                selpl$log2FoldChange) # Indentify candidate miRNA interractively
rownames(selpl)[idx]

rownames(rbind(sh.pl.sb.200[rownames(dfcontr[complete.cases(dfcontr$inheritance200) &
                                               dfcontr$inheritance200 == "PLdominant" & dfcontr$sh.pl.sb.200.log2FoldChange > 1,]),],
               sh.pl.sb.200[rownames(dfcontr[complete.cases(dfcontr$inheritance200) &
                                               dfcontr$inheritance200 == "PLdominant" & dfcontr$sh.pl.sb.200.log2FoldChange < -1,]),]))
write("~/Desktop/Go-analysis/pldominant200.txt",x = rownames(selpl))




########### Gene variability

genetable <- data.frame(gene.id = rownames(mdf),
                        stringsAsFactors = FALSE)

dge <- DGEList(counts = mdf, 
               samples = df.info, 
               genes = genetable)
names(dge)

# Filter gene with less than 10 read count in all samples
keep2 <- filterByExpr(dge, group =dge$samples[,"group"], min.total.count = 10)
dge <- dge[keep2, , keep.lib.sizes=FALSE]

# Normalize
dge <- calcNormFactors(dge)
tmm <- cpm(dge)

write.table(x = tmm, "~/Desktop/mRNA-table-normalized-counts-charr.txt",sep = " ")
write.table(dge$samples[1], "~/Desktop/mRNA-samples-group.txt", sep = " " )

# Compare LCV distributions
lcv<-read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/LCV/gene_noise_github500_all_genes.csv", sep = ",", h = T)
par(mfrow = c(2,2))
hist(lcv$SBxSB_150, main = "a. SBxSB 150ts" , xlab = "LCV")
hist(lcv$PLxPL_150, main = "c. PLxPL 150ts", xlab = "LCV")
hist(lcv$SBxPL_150, main = "b. SBxPL 150ts", xlab = "LCV")
hist(lcv$PLxSB_150, main = "d. PLxSB 150ts", xlab = "LCV")

par(mfrow = c(2,2))
hist(lcv$SBxSB_200, main = "a. SBxSB 200ts", xlab = "LCV")
hist(lcv$PLxPL_200, main = "c. PLxPL 200ts", xlab = "LCV")
hist(lcv$SBxPL_200, main = "b. SBxPL 200ts", xlab = "LCV")
hist(lcv$PLxSB_200, main = "d. PLxSB 200ts", xlab = "LCV")

# Heatmap
lcv.m<-as.matrix(lcv[complete.cases(lcv),2:9])
colnames(lcv.m)<-colnames(lcv[,2:9])
rownames(lcv.m)<-lcv[complete.cases(lcv),1]

df2 <- data.frame(group = as.factor(colnames(lcv.m)))
ann_colors2 = list(
  group =c(PLxPL_150 = "#003300", SBxSB_150 = "#336600", PLxSB_150 = "#99FF99",SBxPL_150 = "#99FFCC",
           PLxPL_200 = "#003366", SBxSB_200 ="#0066CC", PLxSB_200 ="#99CCFF", SBxPL_200 = "#CCE5FF")
)
hmap<- pheatmap(lcv.m)

col_fun = colorRamp2(c(0, 50, 100), c("darkblue", "grey100", "darkgreen"))


ha =rowAnnotation()
ht = Heatmap(lcv.m,km = 10, col = col_fun, name = "LCV", show_row_names = F)
ht2 <- draw(ht)

clist<- list()
for(i in 1:10){
  c4 <- t(t(row.names(lcv.m[row_order(ht2)[[i]],])))
  df.c4 <-data.frame("gene" = t(t(row.names(lcv.m[row_order(ht2)[[i]],]))), lcv.m[row_order(ht2)[[i]],])
  df.c4 <- df.c4 %>% pivot_longer(cols = 2:9 , names_to = "Group", values_to = "LCV", names_ptypes = factor())
  df.c4<-as.data.frame(df.c4)
  df.c4$Group<-as.factor(df.c4$Group)
  df.c4$gene<-as.factor(df.c4$gene)
  clist[[i]]<- df.c4
  write.csv(df.c4$gene, paste("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/LCV/genes-cluster",i,".txt",sep=""))
}
names(clist) <- 1:10

p1<-list(
  R=list(V=1, nu = 0.002))
itr<-1

dfmc <- list()
for(i in 1:10){
  mc4<-MCMCglmm(LCV~Group -1, prior = p1, rcov = ~units, data = clist[[i]],
                nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr)
  dfmc[[i]]<-data.frame(Group = rownames(summary(mc4)[["solutions"]]),
                    pmode= posterior.mode(mc4$Sol),
                    HPDinterval(mc4$Sol)) %>% separate(Group, c(NA,"Group"), sep = 5)
}

for(i in 1:10) {
  write.table(sample(x = clist[[i]]$gene[!duplicated(clist[[i]]$gene)],200),
                                 paste("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/LCV/Results/LCV-genes-location/cluster",i,".txt",sep =""))
  }  
write.csv(as.data.frame(sample(x = clist[[1]]$gene[!duplicated(clist[[2]]$gene)],25))[,2],"~/Documents/.txt",sep ="")
# Plot LCV posterior modes and 95% CrIs

dfmc.p<- data.frame(rbind(dfmc[[1]],dfmc[[2]],dfmc[[3]],dfmc[[4]],dfmc[[5]],
               dfmc[[6]],dfmc[[7]],dfmc[[8]],dfmc[[9]],dfmc[[10]]),
               Cluster = c(rep("Cluster 1",nrow(dfmc[[1]])),rep("Cluster 2",nrow(dfmc[[1]])),
                           rep("Cluster 3",nrow(dfmc[[1]])),rep("Cluster 4",nrow(dfmc[[1]])),
                           rep("Cluster 5",nrow(dfmc[[1]])),rep("Cluster 6",nrow(dfmc[[1]])),
                           rep("Cluster 7",nrow(dfmc[[1]])),rep("Cluster 8",nrow(dfmc[[1]])),
                           rep("Cluster 9",nrow(dfmc[[1]])),rep("Cluster 10",nrow(dfmc[[1]]))))

ggplot(dfmc.p, aes(x=Group, y=pmode)) + 
    geom_point(size = 0.0001) +
    facet_wrap("Cluster",nrow = 5, ncol = 2, scales = "free_y") +
    #ggtitle(paste("Cluster",9)) +
    ylab("LCV") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          title = element_text(size =7))



#### Genes locations on linkage groups

df.lg<-data.frame(rbind(read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/LCV/Results/LCV-genes-location/coordinates/ncbi-cluster1.tsv",h=T),
                 read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/LCV/Results/LCV-genes-location/coordinates/ncbi-cluster2.tsv",h=T),
read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/LCV/Results/LCV-genes-location/coordinates/ncbi-cluster3.tsv",h=T),
read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/LCV/Results/LCV-genes-location/coordinates/ncbi-cluster4.tsv",h=T),
read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/LCV/Results/LCV-genes-location/coordinates/ncbi-cluster5.tsv",h=T),
read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/LCV/Results/LCV-genes-location/coordinates/ncbi-cluster6.tsv",h=T),
read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/LCV/Results/LCV-genes-location/coordinates/ncbi-cluster7.tsv",h=T),
read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/LCV/Results/LCV-genes-location/coordinates/ncbi-cluster8.tsv",h=T),
read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/LCV/Results/LCV-genes-location/coordinates/ncbi-cluster9.tsv",h=T),
read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/LCV/Results/LCV-genes-location/coordinates/ncbi-cluster10.tsv",h=T)),
cluster = c(rep("Cluster1",nrow(dfmc[[1]])),rep("Cluster2",nrow(dfmc[[1]])),
            rep("Cluster3",nrow(dfmc[[1]])),rep("Cluster4",nrow(dfmc[[1]])),
            rep("Cluster5",nrow(dfmc[[1]])),rep("Cluster6",nrow(dfmc[[1]])),
            rep("Cluster7",nrow(dfmc[[1]])),rep("Cluster8",nrow(dfmc[[1]])),
            rep("Cluster9",nrow(dfmc[[1]])),rep("Cluster10",nrow(dfmc[[1]]))))
plot(as.factor(df.lg$cluster),as.factor(df.lg$Chromosome))

