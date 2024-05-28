library(ggplot2)
library(ade4)   # use to PcoA
library(vegan)  # use to calculate the distances
#(1)taxonomy/function-based distance
# Create a dataframe for the points
points_df <- data.frame(x = c(0, 0, 0), y = c(0, 0, 0), Sample = c("S1", "S2", "S3"))
# Create a scatter plot
scatter_plot <- ggplot(data = points_df, aes(x = x, y = y, label = Sample)) +
  geom_point(aes(shape = Sample), size = 6) +
  scale_shape_manual(values = c(0, 1, 2)) +
  scale_color_manual(values = c("orange", "yellow", "blue")) +
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  ggtitle("Taxonomy/Function-based distance")+ coord_fixed()+
  theme_bw()+ xlim(-1, 1) + ylim(-1, 1)
# Display the scatter plot
print(scatter_plot)

#(2)'protein-based distance
data <- read.table("ProteinGroup Table.txt", sep = "\t", header = T, row.names = 1)
distances <- dist(t(data))
distances
pcoa_result =  dudi.pco(distances,
                        scannf = F,   
                        nf=2)         
pcoa_result_data = pcoa_result$li
# Create a dataframe for plotting
plot_df <- data.frame(Sample = rownames(pcoa_result_data), 
                      Axis1 = pcoa_result_data[, 1], 
                      Axis2 = pcoa_result_data[, 2])
# Create an PCoA plot using ggplot2 with different shapes for the three data points
pcoa_plot <- ggplot(data = plot_df, aes(x = Axis1, y = Axis2, label = Sample)) +
  geom_point(aes(shape = Sample), size = 6) +
  scale_shape_manual(values = c(0, 1, 2)) +
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  ggtitle("Protein-based distance")+ coord_fixed()+
  theme_bw()+ xlim(-1.5, 1.5) + ylim(-1.5, 1.5)
# Display the PCoA plot
print(pcoa_plot)

#(3)PhyloFun distance
dist <- read.csv("PhyloFun_distance.csv", row.names=1, header=TRUE)
dist=as.dist(dist)
pcoa_result =  dudi.pco(dist,
                        scannf = F,   
                        nf=2)         
pcoa_result_data = pcoa_result$li
# Create a dataframe for plotting
plot_df <- data.frame(Sample = rownames(pcoa_result_data), 
                      Axis1 = pcoa_result_data[, 1], 
                      Axis2 = pcoa_result_data[, 2])
# Create an PCoA plot using ggplot2 with different shapes for the three data points
pcoa_plot <- ggplot(data = plot_df, aes(x = Axis1, y = Axis2, label = Sample)) +
  geom_point(aes(shape = Sample), size = 6) +
  scale_shape_manual(values = c(0, 1, 2)) +
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  ggtitle("PhyloFunc distance")+ coord_fixed()+
  theme_bw()+ xlim(-0.3, 0.3) + ylim(-0.3, 0.3)
# Display the PCoA plot
print(pcoa_plot)
