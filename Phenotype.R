library(tidyr)
library(ggplot2)
library(MCMCglmm)
library(tidybayes)
library(gridExtra)

d1 <- read.csv("~/RA-phenotype-r1.csv")
d2 <- read.csv("~/RA-phenotype-r2.csv")
d1 <- rbind(d1,d2)

d1 <- d1[, c("Label","Length","Angle")]
d1$Label <- as.factor(d1$Label)

d1$Measure <- rep(c("Articular","Dentary"),nrow(d1)/2)

d1 <- pivot_wider(d1, names_from = "Measure", values_from = c("Length","Angle"))
d1 <- as.data.frame(d1)
d1$Ratio <- d1$Length_Articular/d1$Length_Dentary
d1$Angle <- 360-((d1$Angle_Dentary+180)-d1$Angle_Articular)
d1$Label2 <- d1$Label
d1 <- separate(d1, Label2, into = c(NA, "Cross",NA), c(14,18))


d1 <- separate(d1, Label, into = c(NA, "ind",NA), c(12,24))
d1 <- separate(d1, ind, into = c("ind",NA), "[.]")
d1$fam <- d1$ind 
d1 <- separate(d1, fam, into = c(NA, "fam",NA), c(2,8))

d1$fam.short <- d1$ind 
d1 <- separate(d1, fam.short, into = c(NA, "fam.short",NA), c(6,8))

d3 <- read.table("~/LiekePonsioen_sizes.txt",h=T,numerals = "no.loss") # Dataset of fork length sizes
d3$ind <- paste(d3$Family, formatC(d3$name, width=3, flag="0"), sep="_")
d1<- merge(d1,d3[, c("ind","Length")], by = "ind",all.x = T)
d1$alpha <- 180-(d1$Angle - 90)-90
d1$sin.alpha <- sin(d1$alpha)

############# Articular/Dentary ratio #############

### Boxplot
ggplot(d1, aes(x=Cross, y=Ratio, fill = Cross)) +
  geom_boxplot() + 
  geom_point(position=position_jitter(w=0.12,h=0)) + 
  scale_fill_manual(values = c("darkgreen","lightgreen","cornflowerblue"))+ 
  theme_classic() +
  theme(legend.position = "none",
        title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 14))

### Boxplot
ggplot(d1, aes(x=Length, y=Ratio, color = Cross)) +
  geom_point(position=position_jitter(w=0.12,h=0), size = 0.5 ) + 
  scale_color_manual(values = c("darkgreen","lightgreen","cornflowerblue")) + 
  theme_classic() +
  theme(legend.position = "none",
        title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 14))

############# Articular,Dentary Angle #############

ggplot(d1, aes(x=Cross, y=Angle, fill = Cross)) +
  geom_boxplot() + 
  geom_point(position=position_jitter(w=0.12,h=0)) + 
  scale_fill_manual(values = c("darkgreen","lightgreen","cornflowerblue"))+ 
  theme_classic() +
  theme(legend.position = "none",
        title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 14))

######## Linear Model 

prior1<-list(G=list(
  G1=list(V=1,nu=0.002)),R=list(V=diag(3),n=1))
iter <- 50
m.ratio1<-MCMCglmm(Ratio~Cross - 1,
                           random =~ fam, prior = prior1, 
                           rcov = ~idh(Cross):units, data = d1,verbose = T,
                           nitt = 13000*iter,thin = 10*iter, burnin = 3000*iter)

m.angle1<-MCMCglmm(alpha~Cross - 1,
                   random =~ fam, prior = prior1, 
                   rcov = ~idh(Cross):units, data = d1,verbose = T,
                   nitt = 13000*iter,thin = 10*iter, burnin = 3000*iter)

## Mechanc Advantage

### Fixed effects 

df.ratio <- data.frame(Cross =c(rep("PLPL",nrow(m.ratio1$Sol)),
                                rep("PLSB",nrow(m.ratio1$Sol)),
                                rep("SBSB",nrow(m.ratio1$Sol))),
                       MA = c(m.ratio1$Sol[,"CrossPLPL"],
                              m.ratio1$Sol[,"CrossPLSB"],
                              m.ratio1$Sol[,"CrossSBSB"]))


p.ratio <- ggplot(df.ratio, aes(y = Cross, x = MA)) +  
  geom_halfeyeh(.width = 0.95,fill = "lightblue", color = "#003333",point_interval = mode_hdci) +
  #geom_vline(xintercept = 0.0, color = "lightblue",linetype = "longdash",size = 0.8) +
  ggtitle("(a)") +
  #scale_y_discrete(labels = expression(F[1],PL,SB)) + 
  theme_classic() + 
  xlab("I/O ratio") + ylab("Cross type") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 17, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14,face = "bold"),
        strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text.x = element_text(size = 12,face = "bold", hjust = 0.05))

### Residuals 

df.ratio.r <- data.frame(Cross =c(rep("PLPL",nrow(m.ratio1$VCV)),
                                rep("PLSB",nrow(m.ratio1$VCV)),
                                rep("SBSB",nrow(m.ratio1$VCV))),
                       MA = c(m.ratio1$VCV[,"CrossPLPL.units"],
                              m.ratio1$VCV[,"CrossPLSB.units"],
                              m.ratio1$VCV[,"CrossSBSB.units"]))


p.ratio.r <- ggplot(df.ratio.r, aes(y = Cross, x = MA)) +  
  geom_halfeyeh(.width = 0.95,fill = "lightblue", color = "#003333",point_interval = mode_hdci) +
  #geom_vline(xintercept = 0.0, color = "lightblue",linetype = "longdash",size = 0.8) +
  ggtitle("(b)") +
  #scale_y_discrete(labels = expression(F[1],PL,SB)) + 
  theme_classic() + 
  xlab("Residuals I/O ratio") + ylab("") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 17, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14,face = "bold"),
        strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text.x = element_text(size = 12,face = "bold", hjust = 0.05))

########### Angle 
df.angle <- data.frame(Cross =c(rep("PLPL",nrow(m.angle1$Sol)),
                                rep("PLSB",nrow(m.angle1$Sol)),
                                rep("SBSB",nrow(m.angle1$Sol))),
                       MA = c(m.angle1$Sol[,"CrossPLPL"],
                              m.angle1$Sol[,"CrossPLSB"],
                              m.angle1$Sol[,"CrossSBSB"]))


p.angle <- ggplot(df.angle, aes(y = Cross, x = MA)) +  
  geom_halfeyeh(.width = 0.95,fill = "lightblue", color = "#003333",point_interval = mode_hdci) +
  #geom_vline(xintercept = 0.0, color = "lightblue",linetype = "longdash",size = 0.8) +
  ggtitle("(c)") +
  #scale_y_discrete(labels = expression(F[1],PL,SB)) + 
  theme_classic() + 
  xlab(expression(paste(alpha," ","angle (°)"))) + ylab("Cross type") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 17, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14,face = "bold"),
        strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text.x = element_text(size = 12,face = "bold", hjust = 0.05))

### Residuals 

df.angle.r <- data.frame(Cross =c(rep("PLPL",nrow(m.angle1$VCV)),
                                  rep("PLSB",nrow(m.angle1$VCV)),
                                  rep("SBSB",nrow(m.angle1$VCV))),
                         Angle = c(m.angle1$VCV[,"CrossPLPL.units"],
                                m.angle1$VCV[,"CrossPLSB.units"],
                                m.angle1$VCV[,"CrossSBSB.units"]))


p.angle.r <- ggplot(df.angle.r, aes(y = Cross, x = Angle)) +  
  geom_halfeyeh(.width = 0.95,fill = "lightblue", color = "#003333",point_interval = mode_hdci) +
  #geom_vline(xintercept = 0.0, color = "lightblue",linetype = "longdash",size = 0.8) +
  ggtitle("(d)") +
  #scale_y_discrete(labels = expression(F[1],PL,SB)) + 
  theme_classic() + 
  xlab("Residuals Angle (O,I)") + ylab("") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 17, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14,face = "bold"),
        strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text.x = element_text(size = 12,face = "bold", hjust = 0.05))
grid.arrange(p.ratio,p.ratio.r,p.angle,p.angle.r, ncol= 2)
