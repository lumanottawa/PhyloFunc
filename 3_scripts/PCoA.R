library(vegan)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ade4)    
library(vegan)
library(gridExtra)
# Draw a PcoA plot based on the resulting PhyloFunc distance matrix of mouse gut microbiome
dist <- read.csv("PhyloFunc_distance_mouse.csv", row.names=1, header=TRUE)
dist=as.dist(dist)
pcoa =  dudi.pco(dist,
                 scannf = F,   
                 nf=10)         
data = pcoa$li
data$Fiber_group <- substr(rownames(data), 15, 18)
data$Species_count <- substr(rownames(data), 20, 25)
eig <- pcoa$eig
plot_pif<-ggplot(data=data,aes(x=A1,y=A2,shape=Species_count,color=Fiber_group))+
  geom_point(alpha=1,size=4)+
  scale_color_manual(values = c("HiSF" = "#5F5790",
                                "PEFi" = "#B1764A")) +
  scale_fill_manual(values = c("HiSF" = "#5F5790",
                               "PEFi" = "#B1764A")) +
  stat_ellipse(aes(fill=Fiber_group),type = "norm",geom="polygon",alpha=0.2,color=NA)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  coord_fixed(ratio=1)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=12),axis.title.y = element_text(colour = 'black',size=12),legend.text = element_text(size = 12))+
  ggtitle("PhyloFunc Distance")
ggsave("PhyloFunc_PCoA.pdf", plot = plot_pif, width = 6, height = 4)