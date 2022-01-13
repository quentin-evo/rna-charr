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
write.csv("~/Desktop/miRNA-maternal-test.txt", x = chitable )

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
summary(dfcontr$inheritance150)

# Data frame
dominance <- data.frame(
  type = c(rep("DE150",nrow(sh.pl.sb.150[complete.cases(sh.pl.sb.150$padj) & sh.pl.sb.150$padj < 0.1,])),
           # List of differential expressed genes between pure-morph crosses at 150ts
           rep("DE200",nrow(sh.pl.sb.200[complete.cases(sh.pl.sb.200$padj) & sh.pl.sb.200$padj < 0.1,])),
           # List of differential expressed genes between pure-morph crosses at 200ts
           rep("PLdom150", length(rownames(dfcontr[dfcontr$inheritance150 %in% "PLdominant",]))),
           # List of PL dominant genes at 150ts
           rep("Maternal200", length(rownames(dfcontr[dfcontr$inheritance200 %in% "Maternal",])))
           # List of "maternal" genes at 200ts
           ),
    names = c(rownames(sh.pl.sb.150[complete.cases(sh.pl.sb.150$padj) & sh.pl.sb.150$padj < 0.1,]),
              rownames(sh.pl.sb.200[complete.cases(sh.pl.sb.200$padj) & sh.pl.sb.200$padj < 0.1,]),
              rownames(dfcontr[dfcontr$inheritance150 %in% "PLdominant",]),
              rownames(dfcontr[dfcontr$inheritance200 %in% "Maternal",])))

write.csv(dominance, "~/Desktop/dominant.miRNA.txt")
str(sh.pl.sb.200[sh.pl.sb.200$padj < 0.1 ,])

######## Sort miRNAs 
library(tidyverse)
domi2 <- dominance
domi2$names <- gsub('^....','', domi2$names)
domi2$names <-gsub('*[.][[:digit:]]', '', domi2$names)
domi.t <-table(domi2$names,domi2$type)
domi.t <- as.data.frame(domi.t)
domi.t <- domi.t[domi.t$Freq != 0 ,]

### Only with miRNA families
domi3 <- domi2
domi3$names <-gsub('*-.p', '', domi3$names)
domi3$names <-gsub('*[abcdgj]', '', domi3$names)
domi.t2 <-table(domi3$names,domi3$type)
domi.t2 <- as.data.frame(domi.t2)
domi.t2 <- domi.t2[domi.t2$Freq != 0 ,]
write.csv(domi.t2, "~/Documents/Analyses/Post-zygotic mechanisms/RNAseq/hybridsmiRNAseq/quentin/results/dominant-table.miRNA.txt")
