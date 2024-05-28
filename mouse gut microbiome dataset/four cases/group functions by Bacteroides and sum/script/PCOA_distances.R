library(vegan)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ade4)   # use to PcoA
library(vegan)  # use to calculate the distances
setwd("C:/Users/wangl/Desktop/PhyloFun files/mice dataset/exclude all species of Bacteroides")
summarized_df <- read.csv("mice_data_groupby_taxon_COG.csv", header = T)
bac_rows <- summarized_df[grepl("^Bacteroides", summarized_df$taxon), ]
non_bac_rows <- summarized_df[!grepl("^Bacteroides", summarized_df$taxon), ]
bac_rows$taxon="bacteroides"
summarized_bac <- bac_rows %>%
  group_by(taxon,COG_number)  %>%
  summarize(
    across(starts_with("LFQ"), sum),.groups = "drop"
  )
cell_sum_bac<-rbind(summarized_bac,non_bac_rows)
write.csv(cell_sum_bac, "mice_taxon_COG_sum_bac.csv", quote=F, row.names = F)

#1.bray_curtis
bray_curtis_dist <- vegdist(t(cell_sum_bac[, 3:ncol(cell_sum_bac)]), method = "bray")
bray_curtis_frame <- as.data.frame(as.matrix(bray_curtis_dist))
write.csv(bray_curtis_frame, "Bray_curtis.csv", quote=F, row.names = T)
pcoa =  dudi.pco(bray_curtis_dist,
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
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=15),axis.title.y = element_text(colour = 'black',size=15),legend.text = element_text(size = 10))+
  ggtitle("Bray-Curtis Distance")+
  theme(plot.title = element_text(size = 15))
ggsave("Bray_Curtis_sum_bac.pdf", plot = plot, width = 6, height = 4)

#2.Euclidean
euclidean_dist<-vegdist(t(cell_sum_bac[, 3:ncol(cell_sum_bac)]), method = "euclidean")
euclidean_frame <- as.data.frame(as.matrix(euclidean_dist))
write.csv(euclidean_frame, "Euclidean_sum_bac.csv", quote=F, row.names = T)
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
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=15),axis.title.y = element_text(colour = 'black',size=15),legend.text = element_text(size = 10))+
  ggtitle("Euclidean Distance")+
  theme(plot.title = element_text(size = 15))
ggsave("Euclidean_sum_bac.pdf", plot = plot, width = 6, height = 4)


#3. Jaccard
jaccard_dist<-vegdist(t(cell_sum_bac[, 3:ncol(cell_sum_bac)]), method = "jaccard", binary=TRUE)
jaccard_frame <- as.data.frame(as.matrix(jaccard_dist))
write.csv(jaccard_frame, "BinaryJaccard.csv", quote=F, row.names = T)
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
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=15),axis.title.y = element_text(colour = 'black',size=15),legend.text = element_text(size = 10))+
  ggtitle("Binary Jaccard Distance")+
  theme(plot.title = element_text(size = 15))
ggsave("BinaryJaccard_sum_bac.pdf", plot = plot, width = 6, height = 4)

#4.distance of PiF
dist <- read.csv("PhyloFunc_distance_sum_bac.csv", row.names=1, header=TRUE)
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
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=15),axis.title.y = element_text(colour = 'black',size=15),legend.text = element_text(size = 10))+
  ggtitle("PhyloFunc Distance")+
  theme(plot.title = element_text(size = 15))
ggsave("PhyloFunc Distance_sum_bac.pdf", plot = plot)

