library(vegan)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ade4)   # use to PcoA
library(vegan)  # use to calculate the distances
summarized_df <- read.csv("mice_taxon_COG_except_bac.csv", header = T)

#1.bray_curtis
bray_curtis_dist <- vegdist(t(summarized_df[, 3:ncol(summarized_df)]), method = "bray")
bray_curtis_frame <- as.data.frame(as.matrix(bray_curtis_dist))
write.csv(bray_curtis_frame, "mice_Bray_curtis_distance_excpet_bac.csv", quote=F, row.names = T)
pcoa =  dudi.pco(bray_curtis_dist,
                 scannf = F,  
                 nf=10)         # how many dimensions of PCoA information

data = pcoa$li
data$drug <- substr(rownames(data), 15, 18)
data$concentration <- substr(rownames(data), 20, 25)
eig <- pcoa$eig
plot<-ggplot(data=data,aes(x=A1,y=A2,shape=concentration,color=drug))+
  geom_point(alpha=1,size=4)+
  stat_ellipse(aes(fill=drug),type = "norm",geom="polygon",alpha=0.2,color=NA)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  coord_fixed(ratio=1)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))+
  ggtitle("Bray-Curtis Distance")
ggsave("Bray_Curtis_Distance.pdf", plot = plot, width = 6, height = 4)

#2.Euclidean
euclidean_dist<-vegdist(t(summarized_df[, 3:ncol(summarized_df)]), method = "euclidean")
euclidean_frame <- as.data.frame(as.matrix(euclidean_dist))
write.csv(euclidean_frame, "Euclidean distance.csv", quote=F, row.names = T)
pcoa =  dudi.pco(euclidean_dist,
                 scannf = F,  
                 nf=10)         # how many dimensions of PCoA information

data = pcoa$li
data$drug <- substr(rownames(data), 15, 18)
data$concentration <- substr(rownames(data), 20, 25)
eig <- pcoa$eig
plot<-ggplot(data=data,aes(x=A1,y=A2,shape=concentration,color=drug))+
  geom_point(alpha=1,size=4)+
  stat_ellipse(aes(fill=drug),type = "norm",geom="polygon",alpha=0.2,color=NA)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  coord_fixed(ratio=1)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))+
  ggtitle("Euclidean Distance")
ggsave("Euclidean_distance.pdf", plot = plot, width = 6, height = 4)

#3. Jaccard
jaccard_dist<-vegdist(t(summarized_df[, 3:ncol(summarized_df)]), method = "jaccard", binary=TRUE)
jaccard_frame <- as.data.frame(as.matrix(jaccard_dist))
write.csv(jaccard_frame, "BinaryJaccard.csv", quote=F, row.names = T)
pcoa =  dudi.pco(jaccard_dist,
                 scannf = F,  
                 nf=10)         # how many dimensions of PCoA information

data = pcoa$li
data$drug <- substr(rownames(data), 15, 18)
data$concentration <- substr(rownames(data), 20, 25)
eig <- pcoa$eig
plot<-ggplot(data=data,aes(x=A1,y=A2,shape=concentration,color=drug))+
  geom_point(alpha=1,size=4)+
  stat_ellipse(aes(fill=drug),type = "norm",geom="polygon",alpha=0.2,color=NA)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  coord_fixed(ratio=1)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))+
  ggtitle("Binary Jaccard Distance")
ggsave("Binary Jaccard Distance.pdf", plot = plot, width = 6, height = 4)

#4.distance of PiF
dist <- read.csv("PhyloFunc_distance_exclude_Bac.csv", row.names=1, header=TRUE)
dist=as.dist(dist)
pcoa =  dudi.pco(dist,
                 scannf = F,   
                 nf=10)         
data = pcoa$li
data$drug <- substr(rownames(data), 15, 18)
data$concentration <- substr(rownames(data), 20, 25)
eig <- pcoa$eig
plot<-ggplot(data=data,aes(x=A1,y=A2,shape=concentration,color=drug))+
  geom_point(alpha=1,size=4)+
  stat_ellipse(aes(fill=drug),type = "norm",geom="polygon",alpha=0.2,color=NA)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  coord_fixed(ratio=1)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))+
  ggtitle("PhyloFunc Distance")
ggsave("PhyloFunc Distance.pdf", plot = plot)

