library(vegan)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ade4)   # use to PcoA
library(vegan)  # use to calculate the distances

summarized_df <- read.csv("mice_taxon_COG_only Bcateroides.csv", header = T)

#1.bray_curtis
bray_curtis_dist <- vegdist(t(summarized_df[, 3:ncol(summarized_df)]), method = "bray")
bray_curtis_frame <- as.data.frame(as.matrix(bray_curtis_dist))
write.csv(bray_curtis_frame, "Bray_curtis distance.csv", quote=F, row.names = T)
pcoa =  dudi.pco(bray_curtis_dist,
                 scannf = F,  
                 nf=10)         # how many dimensions of PCoA information

data = pcoa$li
data$Diet <- substr(rownames(data), 15, 18)
data$Community_Members <- substr(rownames(data), 20, 25)
eig <- pcoa$eig
plot<-ggplot(data=data,aes(x=A1,y=A2,shape=Community_Members,color=Diet))+
  geom_point(alpha=1,size=4)+
  stat_ellipse(aes(fill=Diet),type = "norm",geom="polygon",alpha=0.2,color=NA)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  coord_fixed(ratio=1)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))+
  ggtitle("Bray-Curtis Distance")+
  theme(plot.title = element_text(size = 15))
ggsave("Bray_Curtis_Bac.pdf", plot = plot, width = 6, height = 4)

#2.Euclidean
euclidean_dist<-vegdist(t(summarized_df[, 3:ncol(summarized_df)]), method = "euclidean")
euclidean_frame <- as.data.frame(as.matrix(euclidean_dist))
write.csv(euclidean_frame, "cell_euclidean_LFQ.csv", quote=F, row.names = T)
pcoa =  dudi.pco(euclidean_dist,
                 scannf = F,  
                 nf=10)         # how many dimensions of PCoA information

data = pcoa$li
data$Diet <- substr(rownames(data), 15, 18)
data$Community_Members <- substr(rownames(data), 20, 25)
eig <- pcoa$eig
plot<-ggplot(data=data,aes(x=A1,y=A2,shape=Community_Members,color=Diet))+
  geom_point(alpha=1,size=4)+
  stat_ellipse(aes(fill=Diet),type = "norm",geom="polygon",alpha=0.2,color=NA)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  coord_fixed(ratio=1)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))+
  ggtitle("Euclidean Distance")+
  theme(plot.title = element_text(size = 15))
ggsave("Euclidean_Bac.pdf", plot = plot, width = 6, height = 4)


#3. Jaccard
jaccard_dist<-vegdist(t(summarized_df[, 3:ncol(summarized_df)]), method = "jaccard", binary=TRUE)
jaccard_frame <- as.data.frame(as.matrix(jaccard_dist))
write.csv(jaccard_frame, "cell_jaccard_LFQ.csv", quote=F, row.names = T)
pcoa =  dudi.pco(jaccard_dist,
                 scannf = F,  
                 nf=10)         # how many dimensions of PCoA information

data = pcoa$li
data$Diet <- substr(rownames(data), 15, 18)
data$Community_Members <- substr(rownames(data), 20, 25)
eig <- pcoa$eig
plot<-ggplot(data=data,aes(x=A1,y=A2,shape=Community_Members,color=Diet))+
  geom_point(alpha=1,size=4)+
  stat_ellipse(aes(fill=Diet),type = "norm",geom="polygon",alpha=0.2,color=NA)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  #coord_fixed(ratio=eig[2]/eig[1])+
  coord_fixed(ratio=1)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))+
  ggtitle("Binary Jaccard Distance")+
  theme(plot.title = element_text(size = 15))
ggsave("Jaccard_Bac.pdf", plot = plot, width = 6, height = 4)

#4.distance of PiF
dist <- read.csv("PhyloFunc_distance_only Bac.csv", row.names=1, header=TRUE)
dist=as.dist(dist)
pcoa =  dudi.pco(dist,
                 scannf = F,   
                 nf=10)         
data = pcoa$li
data$Diet <- substr(rownames(data), 15, 18)
data$Community_Members <- substr(rownames(data), 20, 25)
eig <- pcoa$eig
plot<-ggplot(data=data,aes(x=A1,y=A2,shape=Community_Members,color=Diet))+
  geom_point(alpha=1,size=4)+
  stat_ellipse(aes(fill=Diet),type = "norm",geom="polygon",alpha=0.2,color=NA)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  coord_fixed(ratio=1)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))+
  ggtitle("PhyloFunc Distance")+
  theme(plot.title = element_text(size = 15))
ggsave("PhyloFunc_Bac.pdf", plot = plot)

