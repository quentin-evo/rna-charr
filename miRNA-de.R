load("/Users/Muscardin/Documents/Analyses/Post-zygotic mechanisms/RNAseq/hybridsmiRNAseq/hybridsmiRNAseq/miRNA-de.RData")
library(tidyr)
library("DESeq2")
library("tximport")
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(limma)
library(edgeR)
library(ComplexHeatmap)
library(circlize)
library(MCMCglmm)

# Get miRNA count matrix

d.raw<-read.table("/Users/Muscardin/Documents/Analyses/Post-zygotic mechanisms/RNAseq/hybridsmiRNAseq/quentin/mirna-counts-for-analyses.txt", h=T)
d.samples<-read.table("/Users/Muscardin/Documents/Analyses/Post-zygotic mechanisms/RNAseq/hybridsmiRNAseq/quentin/config.txt")

d.samples<-separate(d.samples,col = "V1", into = c(NA,"V1"),sep = "d/") %>% separate(col = "V1", into = c("V1",NA),sep = "_S")
d.samples<-d.samples[order(d.samples$V2),]
colnames(d.raw[,5:76])<-d.samples$V1

### List of nonredundant miRNAs (from CD-hit EST output)
uniq.prec <- read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/miRNA-target/nonredundant.miRNA.names.txt")
uniq.prec <- as.data.frame(uniq.prec) %>% separate(col = "V1", into = c(NA, "miRNA"), sep = 1)
uniq.mature <- merge(d.raw[,c(1,3)],uniq.prec, by.x = "precursor", by.y="miRNA", all.x = F )
saveRDS(uniq.mature, "~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/miRNA-target/uniq.mature.RDS")

# dataframe filtered from redundant sequences identified with CD-hit
d.raw.filt <- d.raw[d.raw$precursor %in% uniq.mature$precursor,]

dfde<-as.matrix(d.raw.filt[,5:76])
row.names(dfde)<-d.raw.filt[,1]
rownames(dfde)<-make.unique(rownames(dfde))
colnames(dfde)<-d.samples$V1

# Get sample info 
df.info<-read.csv("/Users/Muscardin/Documents/Analyses/Post-zygotic mechanisms/RNAseq/hybridsmiRNAseq/hybridsmiRNAseq/dataInfo/hybridsCharrDataFrameDagny.csv",h=T, sep = ";")
df.info$Sample<-sub("_..._","-",df.info$Sample)
df.info$Batch<-as.factor(df.info$Batch)
df.info$Timepoint<-as.factor(df.info$Timepoint)
row.names(df.info)<-df.info$Sample
df.info<-df.info[,2:6]
df.info$group<-paste(df.info$Cross,df.info$Timepoint,sep = "_")
df.info<-df.info[colnames(dfde),]


# Data set
dds <- DESeqDataSetFromMatrix(countData = dfde,
                              colData = df.info,
                              design = ~Batch + group)
# miRNAs with less than 10 reads 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


# Dataset without H1512-200B
dds2 <- dds[,colnames(dds) != "H1512-200B"]
# Dataset with H1512-200B nor Batch 3
dds3<-dds2[, dds2$Batch != "3"]
dds3$Batch<-droplevels(dds3$Batch)

##### Exploratory analyses 
vsdo<-vst(dds, blind=FALSE) # with outlier 

vsd <- vst(dds2, blind=FALSE) # Transformation by Variance Stabilizing Transformation
#vsdb<-vst(dds3, blind=FALSE) # without Batch3

# Heatmap
select <- order(rowMeans(counts(dds2,normalized=F)),
                decreasing=TRUE)
# Specify colors
ann_colors = list(
  group =c(PLxPL_150 = "#003300", SBxSB_150 = "#336600", PLxSB_150 = "#99FF99",SBxPL_150 = "#99FFCC",
           PLxPL_200 = "#003366", SBxSB_200 ="#0066CC", PLxSB_200 ="#99CCFF", SBxPL_200 = "#CCE5FF"),
  Batch = c("1" = "#202020", "2" = "#808080", "3" = "#E0E0E0")
)
df <- as.data.frame(colData(dds2)[,c("group","Batch")])
pheatmap(assay(vsdo)[select,], cluster_rows=T, show_rownames=FALSE,
         cluster_cols=T, annotation_col=df ,
         annotation_colors = ann_colors
         )


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


# Without batch effect
vsd4<-vsd
assay(vsd4)<-removeBatchEffect(assay(vsd), batch = vsd$Batch, design = model.matrix(~vsd$group))

pcaData2 <- plotPCA(vsd4, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData2, "percentVar"))
pcaData2$group.1<-factor(pcaData2$group.1, levels = c("PLxPL_150","SBxSB_150","PLxSB_150","SBxPL_150","PLxPL_200","SBxSB_200","PLxSB_200","SBxPL_200"))
ggplot(pcaData2, aes(PC1, PC2, color=group.1)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("#003300", "#336600", "99FF99","#99FFCC",
                                "#003366", "#0066CC", "#99CCFF", "#CCE5FF")) +
  labs(color= "Cross & Age") +
  theme_bw()

##### DE analysis
ds_de<-DESeq(dds2)
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
dds2$group<-relevel(dds2$group,"SBxSB_200")
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

# DE tables 


for (i in 1:length(de.list)) {
 write.table(
    as.data.frame(de.list[[i]][order(de.list[[i]]$padj) & complete.cases(de.list[[i]]$padj) & de.list[[i]]$padj < 0.05,]), 
           file = paste("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/hybridsmiRNAseq/quentin/results/", de.names[i],".txt", sep = ""),
    quote = F, 
    sep = '\t')
}


# Outliers
assays(dds)[["cooks"]]

#### MA plot
# PLxPL vs SBxSB 
par(mfrow = c(1,2))
DESeq2::plotMA(sh.pl.sb.150, ylim=c(-2,2)) 
title(main="a.      SBxSB/PLxPL (150ts)"
      , adj = 0
      )
DESeq2::plotMA(sh.pl.sb.200, ylim=c(-2,2), ylab ="")
title(main="b.      SBxSB/PLxPL (200ts)"
      , adj = 0,
      )

# hybrids 150ts
par(mfrow = c(2,2))
DESeq2::plotMA(sh.plsb.plpl.150, ylim=c(-2,2))
title(main="a.  PLxSB/PLxPL (150ts)", adj = 0)
DESeq2::plotMA(sh.plsb.sbsb.150, ylim=c(-2,2))
title(main="c.  PLxSB/SBxSB (150ts)", adj = 0)
DESeq2::plotMA(sh.sbpl.plpl.150, ylim=c(-2,2))
title(main="b.  SBxPL/PLxPL (150ts)", adj = 0)
DESeq2::plotMA(sh.sbpl.sbsb.150, ylim=c(-2,2))
title(main="d.  SBxPL/SBxSB (150ts)", adj = 0)

# hybrids 200ts
par(mfrow = c(2,2))
DESeq2::plotMA(sh.plsb.plpl.200, ylim=c(-2,2))
title(main="a.  PLxSB/PLxPL (200ts)", adj = 0)
DESeq2::plotMA(sh.plsb.sbsb.200, ylim=c(-2,2))
title(main="c.  PLxSB/SBxSB (200ts)", adj = 0)
DESeq2::plotMA(sh.sbpl.plpl.200, ylim=c(-2,2))
title(main="b.  SBxPL/PLxPL (200ts)", adj = 0)
DESeq2::plotMA(sh.sbpl.sbsb.200, ylim=c(-2,2))
title(main="d.  SBxPL/SBxSB (200ts)", adj = 0)

# PLxSB vs. SBxPL
par(mfrow = c(1,2))
DESeq2::plotMA(sh.plsb.sbpl.150, ylim=c(-2,2)) 
title(main="a.      PLxSB/SBxPL (150ts)"
      , adj = 0
)
DESeq2::plotMA(sh.plsb.sbpl.200, ylim=c(-2,2), ylab ="")
title(main="b.       PLxSB/SBxPL (200ts)"
      , adj = 0,
)

idx <- identify(sh.pl.sb.150$baseMean, sh.pl.sb.150$log2FoldChange) # Indentify candidate miRNA interractively
rownames(sh.pl.sb.150)[idx]

# Check for outliers using Cook's distances
assays(dds)[["Cooks"]]

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)


# ##### DE without Batch 3
# ds_de3<-DESeq(dds3)
# resultsNames(ds_de3)
# # 150ts
# r.pl.sb.150.2 <- results(ds_de3,contrast =  c("group","SBxSB_150","PLxPL_150"))
# r.plsb.plpl.150.2 <- results(ds_de3,contrast =  c("group","PLxSB_150","PLxPL_150"))
# r.sbpl.plpl.150.2 <- results(ds_de3,contrast =  c("group","SBxPL_150","PLxPL_150"))
# r.plsb.sbsb.150.2 <- results(ds_de3,contrast =  c("group","PLxSB_150","SBxSB_150"))
# r.sbpl.sbsb.150.2 <- results(ds_de3,contrast =  c("group","SBxPL_150","SBxSB_150"))
# 
# # 200ts
# r.pl.sb.200.2  <- results(ds_de3,contrast =  c("group","SBxSB_200","PLxPL_200"))
# r.plsb.plpl.200.2  <- results(ds_de3,contrast =  c("group","PLxSB_200","PLxPL_200"))
# r.sbpl.plpl.200.2  <- results(ds_de3,contrast =  c("group","SBxPL_200","PLxPL_200"))
# r.plsb.sbsb.200.2  <- results(ds_de3,contrast =  c("group","PLxSB_200","SBxSB_200"))
# r.sbpl.sbsb.200.2  <- results(ds_de3,contrast =  c("group","SBxPL_200","SBxSB_200"))

# Apply shrinkage
dds3$group<-relevel(dds3$group,"SBxSB_200")
ds.de.2<-nbinomWaldTest(ds_de3)

# 150ts
lfcShrink(ds_de3, type="ashr", res = r.pl.sb.150.2, parallel = T,BPPARAM=MulticoreParam(4))
sh.plsb.plpl.150.2 <- lfcShrink(ds_de3, type="ashr", res = r.plsb.plpl.150.2, parallel = T,BPPARAM=MulticoreParam(4))
sh.sbpl.plpl.150.2  <- lfcShrink(ds_de3, type="ashr", res = r.sbpl.plpl.150.2, parallel = T,BPPARAM=MulticoreParam(4))
sh.plsb.sbsb.150.2 <- lfcShrink(ds_de3, type="ashr", res = r.plsb.sbsb.150.2, parallel = T,BPPARAM=MulticoreParam(4))
sh.sbpl.sbsb.150.2 <- lfcShrink(ds_de3, type="ashr", res = r.sbpl.sbsb.150.2, parallel = T,BPPARAM=MulticoreParam(4))

# 200ts
sh.pl.sb.200.2  <- lfcShrink(ds_de3, type="ashr", res = r.pl.sb.200.2, parallel = T,BPPARAM=MulticoreParam(4))
sh.plsb.plpl.200.2 <- lfcShrink(ds_de3, type="ashr", res = r.plsb.plpl.200.2, parallel = T,BPPARAM=MulticoreParam(4))
sh.sbpl.plpl.200.2  <- lfcShrink(ds_de3, type="ashr", res = r.sbpl.plpl.200.2, parallel = T,BPPARAM=MulticoreParam(4))
sh.plsb.sbsb.200.2 <- lfcShrink(ds_de3, type="ashr", res = r.plsb.sbsb.200.2, parallel = T,BPPARAM=MulticoreParam(4))
sh.sbpl.sbsb.200.2 <- lfcShrink(ds_de3, type="ashr", res = r.sbpl.sbsb.200.2, parallel = T,BPPARAM=MulticoreParam(4))

#DE tables
de.list2<-list(sh.pl.sb.150.2,sh.plsb.plpl.150.2, sh.sbpl.plpl.150.2 ,sh.plsb.sbsb.150.2,
              sh.sbpl.sbsb.150.2, sh.pl.sb.200.2,sh.plsb.plpl.200.2, sh.sbpl.plpl.200.2,
              sh.plsb.sbsb.200.2, sh.sbpl.sbsb.200.2)
de.names<- c("sh-pl-sb-150","sh-plsb-plpl-150", "sh-sbpl-plpl-150" ,"sh-plsb-sbsb-150",
             "sh-sbpl-sbsb-150", "sh-pl-sb-200","sh-plsb-plpl-200", "sh-sbpl-plpl-200",
             "sh-plsb-sbsb-200", "sh-sbpl-sbsb-200")

for (i in 1:length(de.list2)) {
  write.csv(
    as.data.frame(de.list2[[i]][order(de.list2[[i]]$padj) & complete.cases(de.list2[[i]]$padj) & de.list2[[i]]$padj < 0.05,]),
    file = paste("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/hybridsmiRNAseq/quentin/results/without-batch3/", de.names[i],".csv", sep=""))
}


# Outliers
assays(dds)[["cooks"]]

#### MA plot
# PLxPL vs SBxSB 
par(mfrow = c(1,2))
plotMA(sh.pl.sb.150.2, ylim=c(-2,2)) 
title(main="a.      PLxPL/SBxSB (150ts)"
      , adj = 0
)
plotMA(sh.pl.sb.200.2, ylim=c(-2,2), ylab ="")
title(main="b.      PLxPL/SBxSB (200ts)"
      , adj = 0,
)

# hybrids 150ts
par(mfrow = c(2,2))
plotMA(sh.plsb.plpl.150.2, ylim=c(-2,2))
title(main="a.  PLxSB/PLxPL (150ts)", adj = 0)
plotMA(sh.plsb.sbsb.150.2, ylim=c(-2,2))
title(main="c.  PLxSB/SBxSB (150ts)", adj = 0)
plotMA(sh.sbpl.plpl.150.2, ylim=c(-2,2))
title(main="b.  SBxPL/PLxPL (150ts)", adj = 0)
plotMA(sh.sbpl.sbsb.150.2, ylim=c(-2,2))
title(main="d.  SBxPL/SBxSB (150ts)", adj = 0)

# hybrids 200ts
par(mfrow = c(2,2))
DESeq::plotMA(sh.plsb.plpl.200.2, ylim=c(-2,2))
title(main="a.  PLxSB/PLxPL (200ts)", adj = 0)
DESeq::plotMA(sh.plsb.sbsb.200.2, ylim=c(-2,2))
title(main="c.  PLxSB/SBxSB (200ts)", adj = 0)
DESeq::plotMA(sh.sbpl.plpl.200.2, ylim=c(-2,2))
title(main="b.  SBxPL/PLxPL (200ts)", adj = 0)
DESeq::plotMA(sh.sbpl.sbsb.200.2, ylim=c(-2,2))
title(main="d.  SBxPL/SBxSB (200ts)", adj = 0)

idx <- identify(sh.sbpl.sbsb.150.2$baseMean, sh.sbpl.sbsb.150.2$log2FoldChange) # Indentify candidate miRNA interractively
rownames(sh.sbpl.sbsb.150.2)[idx]


# Check for outliers using Cook's distances
assays(dds)[["Cooks"]]

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)



##### Plot gene counts

############### miR-34
p.miR34 <- plotCounts(dds2, gene="dre-miR-34a", intgroup=c("group","Batch"), 
                      transform = F , normalized = T , returnData=TRUE)
p.miR34$log.count<-log(p.miR34$count)
p.miR34<-separate(p.miR34,"group",c("Cross","Age"),sep="_")

batchnames<-c("Batch 1","Batch 2","Batch 3")
names(batchnames)<-c("1","2","3")

ggplot(p.miR34, aes(x=Age, y=log.count,color=Cross,shape = Batch)) + 
  ggtitle("miR-34") +
  geom_point(position=position_jitter(w=0.12,h=0)) +
  #facet_wrap(~Batch,scales = "free", labeller = labeller(Batch = batchnames)) + 
  scale_color_manual(values = c("darkgreen","lightgreen","lightblue","blue"))+ 
  theme_bw()

# log2 corrected for batch effect
p.miR34$corrected<-assay(vsd4)["dre-miR-34a",]

ggplot(p.miR34, aes(x=Age, y=corrected,color=Cross)) + 
  ggtitle("miR-34") + 
  ylab("Transformed counts") +
  xlab(expression(paste("Relative age (", tau ["s"],")"))) +
  geom_point(position=position_jitter(w=0.12,h=0)) +
  #facet_wrap(~Batch,scales = "free", labeller = labeller(Batch = batchnames)) + 
  scale_color_manual(values = c("darkgreen","lightgreen","lightblue","blue"))+ 
  theme_bw()

plotmR34 <- ggplot(p.miR34[p.miR34$Age == 150,], aes(x=Cross, y=corrected, fill = Cross)) +
  geom_boxplot() + 
  ggtitle("(c) dre-miR-34a") + 
  geom_point(position=position_jitter(w=0.12,h=0)) + 
  scale_fill_manual(values = c("darkgreen","lightgreen","lightblue","cornflowerblue"))+ 
  ylab("Transformed counts") +
  theme_bw() +
  theme(legend.position = "none")

# ###### With SE
# sl <- ddply(p.miR34, c("Cross", "Age"), summarise,
#             N    = length(log.count),
#             mean = mean(log.count,na.rm = T),
#             sd   = sd(log.count,na.rm = T),
#             se   = sd / sqrt(N)
# )
# pd <- position_dodge(0.1)
# ggplot(sl, aes(x=Age, y=mean, colour= Cross,group = Cross)) + 
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black",width=.1, position=pd) +
#   geom_line(position=pd) +
#   geom_point(position=pd,size = 2) + theme_bw() + theme(panel.grid.major = element_blank(),
#                                                         panel.grid.minor = element_blank()) +
#   labs(title="miR-34") + ylab("Mean log(count)") +
#   scale_color_manual(values = c("darkgreen","lightgreen","lightblue","blue"))+ 
#   xlab("Relative age")


#miR-100
p.miR100 <- plotCounts(dds2, gene="dre-miR-100-5p", intgroup=c("group","Batch"), 
                       transform = F , normalized = T , returnData=TRUE)
p.miR100$log.count<-log(p.miR100$count)
p.miR100<-separate(p.miR100,"group",c("Cross","Age"),sep="_")

batchnames<-c("Batch 1","Batch 2","Batch 3")
names(batchnames)<-c("1","2","3")
ggplot(p.miR100, aes(x=Age, y=log.count,color=Cross,shape=Batch)) + 
  ggtitle("miR-100") +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  #facet_wrap(~Batch,scales = "free", labeller = labeller(Batch = batchnames)) + 
  theme_bw() + 
  scale_color_manual(values = c("darkgreen","lightgreen","lightblue","blue"))

# log2 corrected for batch effect
p.miR100$corrected<-assay(vsd4)["dre-miR-100-5p",]

ggplot(p.miR100, aes(x=Age, y=corrected,color=Cross)) + 
  ggtitle("miR-100") +
  geom_point(position=position_jitter(w=0.12,h=0)) +
  #facet_wrap(~Batch,scales = "free", labeller = labeller(Batch = batchnames)) + 
  scale_color_manual(values = c("darkgreen","lightgreen","lightblue","blue"))+ 
  ylab("Transformed counts") +
  theme_bw()

plotmR100 <- ggplot(p.miR100[p.miR100$Age == 150,], aes(x=Cross, y=corrected, fill = Cross)) +
  geom_boxplot() + 
  ggtitle("(a) dre-miR-100-5p") + 
  geom_point(position=position_jitter(w=0.12,h=0)) + 
  scale_fill_manual(values = c("darkgreen","lightgreen","lightblue","cornflowerblue"))+ 
  ylab("Transformed counts") + xlab("") +
  theme_bw() +
  theme(legend.position = "none")

#ssa-miR-8160-5p
p.miR8160 <- plotCounts(dds2, gene="ssa-miR-8160-5p", intgroup=c("group","Batch"), 
                       transform = F , normalized = T , returnData=TRUE)
p.miR8160$log.count<-log(p.miR8160$count)
p.miR8160<-separate(p.miR8160,"group",c("Cross","Age"),sep="_")

batchnames<-c("Batch 1","Batch 2","Batch 3")
names(batchnames)<-c("1","2","3")
ggplot(p.miR8160, aes(x=Age, y=log.count,color=Cross,shape=Batch)) + 
  ggtitle("miR-8160-5p") +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  #facet_wrap(~Batch,scales = "free", labeller = labeller(Batch = batchnames)) + 
  theme_bw() + 
  scale_color_manual(values = c("darkgreen","lightgreen","lightblue","blue"))

# log2 corrected for batch effect
p.miR8160$corrected<-assay(vsd4)["ssa-miR-8160-5p",]

ggplot(p.miR8160, aes(x=Age, y=corrected,color=Cross)) + 
  ggtitle("miR-8160-5p") +
  geom_point(position=position_jitter(w=0.12,h=0)) +
  #facet_wrap(~Batch,scales = "free", labeller = labeller(Batch = batchnames)) + 
  scale_color_manual(values = c("darkgreen","lightgreen","lightblue","blue"))+ 
  ylab("Transformed counts") +
  theme_bw()

plotmR8160 <- ggplot(p.miR8160[p.miR8160$Age == 150,], aes(x=Cross, y=corrected, fill = Cross)) +
  geom_boxplot() + 
  ggtitle("(d) ssa-miR-8160-5p") + 
  geom_point(position=position_jitter(w=0.12,h=0)) + 
  scale_fill_manual(values = c("darkgreen","lightgreen","lightblue","cornflowerblue"))+ 
  ylab("") +
  theme_bw() +
  theme(legend.position = "none")


#ssa-miR-181a-5p
p.miR181 <- plotCounts(dds2, gene="ssa-miR-181a-5p", intgroup=c("group","Batch"), 
                        transform = F , normalized = T , returnData=TRUE)
p.miR181$log.count<-log(p.miR181$count)
p.miR181<-separate(p.miR181,"group",c("Cross","Age"),sep="_")

batchnames<-c("Batch 1","Batch 2","Batch 3")
names(batchnames)<-c("1","2","3")
ggplot(p.miR181, aes(x=Age, y=log.count,color=Cross,shape=Batch)) + 
  ggtitle("ssa-miR-181a-5p") +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  #facet_wrap(~Batch,scales = "free", labeller = labeller(Batch = batchnames)) + 
  theme_bw() + 
  scale_color_manual(values = c("darkgreen","lightgreen","lightblue","blue"))

# log2 corrected for batch effect
p.miR181$corrected<-assay(vsd4)["ssa-miR-181a-5p",]

ggplot(p.miR181, aes(x=Age, y=corrected,color=Cross)) + 
  ggtitle("ssa-miR-181a-5p") +
  geom_point(position=position_jitter(w=0.12,h=0)) +
  #facet_wrap(~Batch,scales = "free", labeller = labeller(Batch = batchnames)) + 
  scale_color_manual(values = c("darkgreen","lightgreen","lightblue","blue"))+ 
  ylab("Transformed counts") +
  theme_bw()

plotmR181 <- ggplot(p.miR181[p.miR181$Age == 150,], aes(x=Cross, y=corrected, fill = Cross)) +
  geom_boxplot() + 
  ggtitle("(b) ssa-miR-181a-5p") + 
  geom_point(position=position_jitter(w=0.12,h=0)) + 
  scale_fill_manual(values = c("darkgreen","lightgreen","lightblue","cornflowerblue"))+ 
  ylab("") + xlab("") +
  theme_bw() +
  theme(legend.position = "none")


grid.arrange(plotmR100, plotmR181,plotmR34, plotmR8160)


######## Gene variability 

rownames(dfde)<-make.unique(rownames(dfde))

genetable <- data.frame(gene.id = rownames(dfde),
                        stringsAsFactors = FALSE)
stopifnot(all(rownames(dfde) == rownames(dds2)))
dge <- DGEList(counts = dfde, 
               samples = df.info, 
               genes = genetable)
names(dge)

# Filter gene with less than 10 read count in all samples
keep2 <- filterByExpr(dge, group =dge$samples[,"group"], min.total.count = 10)
dge <- dge[keep2, , keep.lib.sizes=FALSE]

# Normalize
dge <- calcNormFactors(dge)
tmm <- cpm(dge)
tmm <- as.data.frame(tmm)

# Select non species redundant miRNA
tmm$gene <- rownames(tmm)
uniq.mature <- readRDS("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/miRNA-target/uniq.mature.RDS")
tmm <- merge(tmm, as.data.frame(uniq.mature[,2]), by.x = "gene", by.y = "uniq.mature[, 2]", all = F)
row.names(tmm) <- make.unique(tmm$gene)
tmm <- tmm[, c(2:72)]
write.table(x = tmm, "~/Desktop/table-normalized-counts-charr.txt",sep = " ")
write.table(dge$samples[1], "~/Desktop/samples-group.txt", sep = " " )

# LCV

lcv<-read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/LCV-master/gene_noise_github500_all_genes.csv", sep = ",", h = T)
par(mfrow = c(2,2))
hist(lcv$SBxSB_150, main = "a. SBxSB 150ts" , xlab = "LCV", col = "lightblue")
hist(lcv$PLxPL_150, main = "c. PLxPL 150ts", xlab = "LCV", col = "lightblue")
hist(lcv$SBxPL_150, main = "b. SBxPL 150ts", xlab = "LCV", col = "lightblue")
hist(lcv$PLxSB_150, main = "d. PLxSB 150ts", xlab = "LCV", col = "lightblue")


par(mfrow = c(2,2))
hist(lcv$SBxSB_200, main = "a. SBxSB 200ts", xlab = "LCV", col = "lightblue")
hist(lcv$PLxPL_200, main = "c. PLxPL 200ts", xlab = "LCV", col = "lightblue")
hist(lcv$SBxPL_200, main = "b. SBxPL 200ts", xlab = "LCV", col = "lightblue")
hist(lcv$PLxSB_200, main = "d. PLxSB 200ts", xlab = "LCV", col = "lightblue")
### GGplot histogram in lov-hist.R


lcv.m<-as.matrix(lcv[complete.cases(lcv),2:9])
colnames(lcv.m)<-colnames(lcv[,2:9])
rownames(lcv.m)<-lcv[complete.cases(lcv),1]

df2 <- data.frame(group = as.factor(colnames(lcv.m)))
ann_colors2 = list(
  group =c(PLxPL_150 = "#003300", SBxSB_150 = "#336600", PLxSB_150 = "#99FF99",SBxPL_150 = "#99FFCC",
           PLxPL_200 = "#003366", SBxSB_200 ="#0066CC", PLxSB_200 ="#99CCFF", SBxPL_200 = "#CCE5FF")
)
hmap<- pheatmap(lcv.m)

library("ComplexHeatmap") ## For heatmap


## Creating heatmap with three clusters (See the ComplexHeatmap documentation for more options)

col_fun = colorRamp2(c(0, 50, 100), c("darkblue", "grey100", "darkgreen"))
col_fun(seq(-3, 3))

ha =rowAnnotation()
ht = Heatmap(lcv.m,km = 10, col = col_fun, name = "LCV", show_row_names = F)
ht2 <- draw(ht)
c4 <- t(t(row.names(lcv.m[row_order(ht2)[[4]],])))
View(c4)

# List of non species redundant RNAs
c4.2<-separate(as.data.frame(c4),col = "V1", sep = 4,into = c(NA, "V1"))
c4.2<- c4.2[!duplicated(c4.2$V1),]
length(c4.2)

df.c4 <-data.frame(mature = t(t(row.names(lcv.m[row_order(ht2)[[4]],]))), lcv.m[row_order(ht2)[[4]],]) %>% separate(,col = "mature", sep = 4,into = c(NA, "mature"))
df.c4 <- df.c4[!duplicated(df.c4$mature),] %>% pivot_longer(cols = 2:9 , names_to = "Group", values_to = "LCV", names_ptypes = factor())
df.c4<-as.data.frame(df.c4)
df.c4$Group<-as.factor(df.c4$Group)
df.c4$mature<-as.factor(df.c4$mature)

write.csv(df.c4, "/Users/Muscardin/Documents/Analyses/Post-zygotic mechanisms/RNAseq/hybridsmiRNAseq/quentin/results/LCV/cluster4.csv")


#### Group differences in c4 with MCMCglmm
p1<-list(
  R=list(V=1, nu = 0.002))
itr<-1
mc4<-MCMCglmm(LCV~Group -1, prior = p1, rcov = ~units, data = df.c4,
              nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr)

dfmc4<-data.frame(Group = rownames(summary(mc4)[["solutions"]]),
                  pmode= posterior.mode(mc4$Sol),
                  HPDinterval(mc4$Sol)) %>% separate(Group, c(NA,"Group"), sep = 5)
ggplot(dfmc4, aes(x=Group, y=pmode)) + 
  geom_point() +
  ggtitle("Cluster 4") + 
  ylab("LCV") +
  xlab("Group") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))
 

# Dataframe of clusters and associated miRNAs
clist<- list()
for(i in 1:10){
  c4 <- t(t(row.names(lcv.m[row_order(ht2)[[i]],])))
  df.c4 <-data.frame("gene" = t(t(row.names(lcv.m[row_order(ht2)[[i]],]))), lcv.m[row_order(ht2)[[i]],])
  df.c4 <- df.c4 %>% pivot_longer(cols = 2:9 , names_to = "Group", values_to = "LCV", names_ptypes = factor())
  df.c4<-as.data.frame(df.c4)
  df.c4$Group<-as.factor(df.c4$Group)
  df.c4$gene<-as.factor(df.c4$gene)
  clist[[i]]<- df.c4
}

# Make class for miRNA across species and without 3p - 5p information. 
for(i in 1:10) {
  c1 <- separate(clist[[i]], col = "gene", into = c(NA, "fam"), sep = 4, remove = F)
  c2 <- separate(c1, col = "fam", into = c("fam", NA), sep = "-3p", remove = F)
  clist[[i]] <- c2
}
names(clist) <- 1:10
dfcl.mi <- dplyr::bind_rows(clist, .id = "miRNACluster")
#dfcl.mi <- dfcl.mi[!duplicated(dfcl.mi[,1:2]), 1:2]
saveRDS(dfcl.mi, file = "~/Desktop/dfcl.mi.rds")


### 
p1<-list(G=list(
  G1=list(V=1,nu=0.002)),
  R=list(V=1, nu = 0.002))
itr<-10

dfmc <- list()
for(i in 1:10){
  mc4<-MCMCglmm(LCV~Group -1, random = ~fam, prior = p1, rcov = ~units, data = clist[[i]],
                nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr)
  dfmc[[i]]<-data.frame(Group = rownames(summary(mc4)[["solutions"]]),
                        pmode= posterior.mode(mc4$Sol),
                        HPDinterval(mc4$Sol)) %>% separate(Group, c(NA,"Group"), sep = 5)
}

for(i in 1:10) {
  write.table(clist[[i]]$gene[!duplicated(clist[[i]]$gene)],
              paste("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/",i,".txt",sep =""))
}

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




#### DE & Variability 
cv.de <- list()
for(i in 1:length(names(contr))) {
    cv.de[[i]]<- merge(cbind(contr[[i]][, c("log2FoldChange","padj")],
                             "X" = rownames(contr[[i]]),
                             "contrast" = rep(names(contr)[[i]], nrow(contr[[i]])) ), 
                       lcv, by = "X")
 }
cv.de2<-as.data.frame(rbind(cv.de[[1]],cv.de[[2]],cv.de[[3]],cv.de[[4]],cv.de[[5]],
                            cv.de[[6]],cv.de[[7]],
              cv.de[[8]],cv.de[[9]],cv.de[[10]],cv.de[[11]],cv.de[[12]])) %>% 
  pivot_longer(cols = ends_with("0"), names_to = "cross", values_to = "LCV") %>% 
  separate("cross",into = c("cross","age"),sep = "_") %>% 
  separate("contrast",into = c("contrast","contrast_age"),sep = -3)
cv.de2 <- as.data.frame(cv.de2[cv.de2$contrast_age == cv.de2$age ,])
cv.de2$p.val <-ifelse(cv.de2$padj < 0.05, "<0.05","NS")
cv.de2$abslog2FoldChange  <- abs(cv.de2$log2FoldChange)
cv.de2<- na.omit(cv.de2)


rm(cv.de)

### Remove redundant species
cv.de2<-separate(cv.de2,col = "X", sep = 4,into = c(NA, "V1"))
cv.de2<- cv.de2[!duplicated(cv.de2[, c("V1","cross","age","contrast")]),]

# At 150ts
ggplot(cv.de2[cv.de2$age == 150,], 
       aes(x = LCV, y = abslog2FoldChange, colour = contrast)) + 
  geom_point(aes(shape = p.val)) + scale_shape_manual(values=c(19,3)) +
  facet_wrap(.~cross,nrow = 2,ncol = 2) + 
  labs(title = "150ts",
       y=expression("| "*log[2]*("Fold change")*" |")) +
  theme_bw()

# At 200ts
ggplot(cv.de2[cv.de2$age == 200,], 
       aes(x = LCV, y = abslog2FoldChange, colour = contrast)) + 
  geom_point(aes(shape = p.val)) + scale_shape_manual(values=c(19,3)) +
  facet_wrap(.~cross,nrow = 2,ncol = 2) + 
  labs(title = "200ts",
       y=expression("| "*log[2]*("Fold change")*" |")) +
  theme_bw()



##### export sh.pl.sb.200 for target analyses

saveRDS(data.frame(sh.pl.sb.200[complete.cases(sh.pl.sb.200$padj) ,
                                    c("log2FoldChange","padj")],
           miRNA.names = rownames(sh.pl.sb.200[complete.cases(sh.pl.sb.200$padj) ,])),
        "~/Desktop/miRNA-foldchange.rds")
