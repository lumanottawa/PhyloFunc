# calculate the other 3 distances matrix
#read original data files in the folder \1_database_search_results\3_human gut microbiome:
library(DESeq2)
library(dplyr)
library(vegan)  # use to calculate the distances

data <- read.csv("Protein_LCA_species_COG_sum_V49_LFQ.csv", header = T)
data_dele_12col=data[,3:ncol(data)]
#calculate 3 distances
#1.bray_curtis
bray_curtis_dist <- vegdist(t(data_dele_12col), method = "bray")
bray_curtis_dist <- as.data.frame(as.matrix(bray_curtis_dist))
write.csv(bray_curtis_dist, "Bray_Curtis.csv", quote=F, row.names = T)
#2.Euclidean
euclidean_dist<-vegdist(t(data_dele_12col), method = "euclidean")
euclidean_dist <- as.data.frame(as.matrix(euclidean_dist))
write.csv(euclidean_dist, "Euclidean distance.csv", quote=F, row.names = T)
#3. Jaccard
jaccard_dist<-vegdist(t(data_dele_12col), method = "jaccard", binary=TRUE)
jaccard_dist <- as.data.frame(as.matrix(jaccard_dist))
write.csv(jaccard_dist, "Jaccard distance.csv", quote=F, row.names = T)


