library(vegan)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ade4)   # use to PcoA
library(vegan)  # use to calculate the distances
data <- read.csv("tidy_data.csv", header = T)
#group_by COG_based function data
summarized_df <- data %>%
  group_by(taxon,COG_number) %>%
  summarize(
    across(starts_with("LFQ"), sum)
  )
write.csv(summarized_df, "preprocessed_data.csv", quote=F, row.names = F)
#1.bray_curtis
bray_curtis_dist <- vegdist(t(summarized_df[, 3:ncol(summarized_df)]), method = "bray")
bray_curtis_frame <- as.data.frame(as.matrix(bray_curtis_dist))
write.csv(bray_curtis_frame, "Bray_curtis.csv", quote=F, row.names = T)
pcoa =  dudi.pco(bray_curtis_dist,
                 scannf = F,  
                 nf=10)         

data = pcoa$li
data$diet <- substr(rownames(data), 15, 18)
data$community_members <- substr(rownames(data), 20, 25)
eig <- pcoa$eig
plot<-ggplot(data=data,aes(x=A1,y=A2,shape=community_members,color=diet))+
  geom_point(alpha=1,size=4)+
  stat_ellipse(aes(fill=diet),type = "norm",geom="polygon",alpha=0.2,color=NA)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  coord_fixed(ratio=1)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))+
  ggtitle("Bray_Curtis_distance")
ggsave("Bray_Curtis.pdf", plot = plot, width = 6, height = 4)

#2.Euclidean
euclidean_dist<-vegdist(t(summarized_df[, 3:ncol(summarized_df)]), method = "euclidean")
euclidean_frame <- as.data.frame(as.matrix(euclidean_dist))
write.csv(euclidean_frame, "Euclidean.csv", quote=F, row.names = T)
pcoa =  dudi.pco(euclidean_dist,
                 scannf = F,  
                 nf=10)         # how many dimensions of PCoA information

data = pcoa$li
data$diet <- substr(rownames(data), 15, 18)
data$community_members <- substr(rownames(data), 20, 25)
eig <- pcoa$eig
plot<-ggplot(data=data,aes(x=A1,y=A2,shape=community_members,color=diet))+
  geom_point(alpha=1,size=4)+
  stat_ellipse(aes(fill=diet),type = "norm",geom="polygon",alpha=0.2,color=NA)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  coord_fixed(ratio=1)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))+
  ggtitle("Euclidean_distance")
ggsave("Euclidean distance.pdf", plot = plot, width = 6, height = 4)

#3. Jaccard
jaccard_dist<-vegdist(t(summarized_df[, 3:ncol(summarized_df)]), method = "jaccard", binary=TRUE)
jaccard_frame <- as.data.frame(as.matrix(jaccard_dist))
#write.csv(jaccard_frame, "BinaryJaccard.csv", quote=F, row.names = T)
pcoa =  dudi.pco(jaccard_dist,
                 scannf = F,  
                 nf=10)         # how many dimensions of PCoA information

data = pcoa$li
data$diet <- substr(rownames(data), 15, 18)
data$community_members <- substr(rownames(data), 20, 25)
eig <- pcoa$eig
plot<-ggplot(data=data,aes(x=A1,y=A2,shape=community_members,color=diet))+
  geom_point(alpha=1,size=4)+
  stat_ellipse(aes(fill=diet),type = "norm",geom="polygon",alpha=0.2,color=NA)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  coord_fixed(ratio=1)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))+
  ggtitle("Binary Jaccard distance")
ggsave("BinaryJaccard.pdf", plot = plot, width = 6, height = 4)

#4.distance of PiF
dist <- read.csv("PhyloFun_distance_mice.csv", row.names=1, header=TRUE)
dist=as.dist(dist)
pcoa =  dudi.pco(dist,
                 scannf = F,   
                 nf=10)         
data = pcoa$li
data$diet <- substr(rownames(data), 15, 18)
data$community_members <- substr(rownames(data), 20, 25)
eig <- pcoa$eig
plot<-ggplot(data=data,aes(x=A1,y=A2,shape=community_members,color=diet))+
  geom_point(alpha=1,size=4)+
  stat_ellipse(aes(fill=diet),type = "norm",geom="polygon",alpha=0.2,color=NA)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  coord_fixed(ratio=1)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))+
  ggtitle("PhyloFunc_distance")
ggsave("PhyloFun_distance.pdf", plot = plot)

