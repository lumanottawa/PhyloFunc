#calculate other 3 distance matrix, the other three cases also use the same code to calculate the other 3 distances matrix
library(DESeq2)
library(dplyr)
library(vegan)  # use to calculate the distances
#read original data files in the folder \1_database_search_results\2_mouse gut microbiome:
data <- read.csv("preprocessed_data.csv", header = T)
#group_by COG_based function data
summarized_df <- data %>%
  group_by(taxon,COG_number) %>%
  summarize(
    across(starts_with("LFQ"), sum)
  )
#1.Bray_curtis distance
bray_curtis_dist <- vegdist(t(summarized_df[, 3:ncol(summarized_df)]), method = "bray")
bray_curtis_frame <- as.data.frame(as.matrix(bray_curtis_dist))
write.csv(bray_curtis_frame, "Bray_curtis.csv", quote=F, row.names = T)

#2.Euclidean
euclidean_dist<-vegdist(t(summarized_df[, 3:ncol(summarized_df)]), method = "euclidean")
euclidean_frame <- as.data.frame(as.matrix(euclidean_dist))
write.csv(euclidean_frame, "Euclidean.csv", quote=F, row.names = T)

#3. Jaccard distance
jaccard_dist<-vegdist(t(summarized_df[, 3:ncol(summarized_df)]), method = "jaccard", binary=TRUE)
jaccard_frame <- as.data.frame(as.matrix(jaccard_dist))
write.csv(jaccard_frame, "BinaryJaccard.csv", quote=F, row.names = T)
