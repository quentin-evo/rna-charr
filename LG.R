library(GenomicRanges)
library(GenomicFeatures)
library(chromoMap)


##### Annotation object
tx.genome <- makeTxDbFromGFF("~/Desktop/GCF_002910315.2_ASM291031v2_genomic.gff", format = "gff")
key1 <-rownames(dds)[1:10]

columns(tx.genome)

# De genes
DE.lg <- AnnotationDbi::select(tx.genome, keys = rownames(na.omit(sh.pl.sb.200)[na.omit(sh.pl.sb.200$padj) < 0.1,]), columns= c("TXCHROM","TXSTART", "TXEND"), 
                      keytype="GENEID")
# Sequences alias
seqalias <- read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/DE-analyses/LG/GCF_002910315.2.chromAlias.txt",h =T)
DE.lg <- DE.lg[DE.lg$TXCHROM %in% seqalias$alias, ]
DE.lg <- merge(DE.lg, seqalias[, c(1,4)], by.x = "TXCHROM", by.y = "alias")
DE.lg <- data.frame(DE.lg$GENEID, DE.lg$sequenceName, DE.lg$TXSTART, DE.lg$TXEND)
DE.lg <- DE.lg[!duplicated(DE.lg[1]),]
colnames(DE.lg) <- NULL
DE.lg[,5] <- rep("Differentialy-expressed", nrow(DE.lg))
DE.lg[DE.lg[,1] %in% rownames(dfcontr[dfcontr$inheritance200  %in% "PLdominant",]),5] <- "PL-dominant"
DE.lg[DE.lg[,1] %in% rownames(dfcontr[dfcontr$inheritance200  %in% "Maternal",]),5] <- "Maternal"

###### Chromosome object 
chrom <- read.table("~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/mRNA/DE-analyses/LG/GCF_002910315.2.chrom.sizes.txt",h =T)
chrom <- data.frame(chrom$chrLG4q_1_29, rep(1,nrow(chrom)),chrom$X90519428)
chrom <- chrom[chrom$chrom.chrLG4q_1_29 %in% DE.lg[,2],]
colnames(chrom) <- NULL
DE.lg <- DE.lg[DE.lg[,2] %in% chrom[,1],]
### Chr graph
chromoMap(list(chrom),list(DE.lg),
          data_based_color_map = T,
          data_type = "categorical", 
          data_colors = list(c("cornflowerblue","orange","black")),
          legend = T,
          chr_color = c("lightblue"),
          segment_annotation = F,
          canvas_width = 1300,
          canvas_height = 700,
          chr_width = 5,
          chr_length = 10,
          links.lg_x = 500,
          links.lg_y = 10000,
          text_font_size = c(15,15),
          grid_text_size = 15,
          left_margin = 100,
          label_font = 15
          )

rm(tx.genome)

############# mRNA Clusters 
lg.list <- list()
for(i in 1:length(clist)){
  cl.lg <- AnnotationDbi::select(tx.genome, keys = as.vector(clist[[i]]$gene[!duplicated(clist[[1]]$gene)]), columns= c("TXCHROM","TXSTART", "TXEND"), 
                                 keytype="GENEID")
  cl.lg <- cl.lg[cl.lg$TXCHROM %in% seqalias$alias, ]
  cl.lg <- merge(cl.lg, seqalias[, c(1,4)], by.x = "TXCHROM", by.y = "alias")
  cl.lg <- data.frame(cl.lg$GENEID, cl.lg$sequenceName, cl.lg$TXSTART, cl.lg$TXEND)
  cl.lg <- cl.lg[!duplicated(cl.lg[1]),]
  lg.list[[i]] <- cl.lg
  lg.list[[i]]$cluster <- rep(paste("Cluster",names(clist[i]),sep =""),nrow(cl.lg))
}
lg.df <- do.call(rbind,lg.list)
colnames(lg.df) <- NULL
lg.df <- lg.df[lg.df[,2] %in% chrom[,1],]
chrom <- chrom[chrom[,1] %in% lg.df[,2],]

### Chr graph
colnames(lg.df) <- c("V5","NA","NA.1","NA.2","NA.3")
chromoMap(list(chrom),list(lg.df),
          data_based_color_map = T,
          data_type = "categorical", 
          data_colors = list(c("cornflowerblue","orange","black",
                               "lightgreen","green","darkgreen",
                               "lightblue","blue","darkblue",
                               "yellow")),
          legend = T,
          chr_color = c("lightgrey"),
          segment_annotation = F,
          canvas_width = 1300,
          canvas_height = 700,
          chr_width = 5,
          chr_length = 10,
          links.lg_x = 500,
          links.lg_y = 10000,
          text_font_size = c(15,15),
          grid_text_size = 15,
          left_margin = 100,
          label_font = 15
)



##### Distance to nearest 
neiblist <- matrix(NA,ncol = 6,nrow = length(clist))
colnames(neiblist) <- c("mean-obs","lower-obs","upper-obs","mean-rand","lower-rand","upper-rand")
for(i in 1:length(clist)){
  tr1 <- clist[[i]]$gene[!duplicated(clist[[i]]$gene)]
  tr1 <- tr1[tr1 %in% genes(tx.genome)$gene_id]
  cl1 <- distanceToNearest(genes(tx.genome)[tr1])
  cl1 <- cl1@elementMetadata@listData[["distance"]]
  #hist(log(cl1@elementMetadata@listData[["distance"]]))
  neiblist[i, 1] <- mean(cl1 )
  neiblist[i, 2] <- mean(cl1 ) - (sd(cl1)/sqrt(length(cl1)))
  neiblist[i, 3] <- mean(cl1 ) + (sd(cl1)/sqrt(length(cl1)))
  ## random
  ltemp <-do.call(rbind,clist[-i])
  ltemp <- sample(ltemp$gene[!duplicated(ltemp$gene)],length(clist[[i]]$gene[!duplicated(clist[[i]]$gene)]))
  dtemp <- distanceToNearest(genes(tx.genome)[ltemp])
  ltemp<- dtemp@elementMetadata@listData[["distance"]]
  neiblist[i, 4] <- mean(ltemp)
  neiblist[i, 5] <- mean(ltemp) - (sd(ltemp)/sqrt(length(ltemp)))
  neiblist[i, 6] <- mean(ltemp) + (sd(ltemp)/sqrt(length(ltemp)))
  #hist(log(dtemp@elementMetadata@listData[["distance"]]))
}
neiblist <- rbind(neiblist[,1:3],neiblist[,4:6])
neiblist$cat <- c(rep("observed",length(clist)),rep("random",length(clist)))
