# Import libraries
library(cluster)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(corrplot))

# Create functions to calculate the ideal number of clusters
ideal_number_of_clusters <- function(dissimilarity_matrix, hc_method, method) {
  fviz_nbclust(dissimilarity_matrix, 
               FUN = hcut, 
               hc_method = hc_method, 
               method = method, 
               k.max = 10) +
    labs(subtitle = paste("The", method, "Method"))
}

# Read the data
data <- read.table("data/ABC.blastp", 
                   sep = "\t", 
                   header = FALSE, 
                   comment.char = "#")

# Assign names to the columns
colnames(data) <- c("query", 
                    "subject", 
                    "identity", 
                    "alignment_length", 
                    "mismatches", 
                    "gap_opens", 
                    "q_start", 
                    "q_end", 
                    "s_start", 
                    "s_end", 
                    "evalue", 
                    "bit_score")

# Calculate the normalized similarity
similarity <- select(data, query, subject, bit_score)
similarity <- mutate(similarity, 
                     normalized_bit_score = bit_score / max(data$bit_score))

# Regularize diagonal in matrix
for(row in 1:nrow(similarity)) {
  if(similarity[row, "query"] == similarity[row, "subject"]) {
    similarity[row, "normalized_bit_score"] <- 1
  }
}

# Calculate dissimilarity
dissimilarity <- mutate(similarity, dissimilarity = 1 - normalized_bit_score)

# Create a dissimilarity matrix
dissimilarity_matrix <- dissimilarity %>%
  select(query, subject, dissimilarity) %>%
  spread(key = subject, value = dissimilarity) %>%
  column_to_rownames(var = "query")


# Single
## Calculate the ideal number of clusters
ideal_number_of_clusters(dissimilarity_matrix, "single", "wss")
ideal_number_of_clusters(dissimilarity_matrix, "single", "silhouette")
##ideal_number_of_clusters(dissimilarity_matrix, "single", "gap_stat")

## Perform hierarchical clustering
csin <- hclust(dist(dissimilarity_matrix), method = "single")
plot(csin, hang = -1, main = "Single Dendogram", cex = 0.4)

### Cut the dendrogram such that 3 clusters are produced
rect.hclust(csin, k=5, border=2:4)


# Average
ideal_number_of_clusters(dissimilarity_matrix, "average", "wss")
ideal_number_of_clusters(dissimilarity_matrix, "average", "silhouette")
##ideal_number_of_clusters(dissimilarity_matrix, "average", "gap_stat")

## Perform hierarchical clustering
cave <- hclust(dist(dissimilarity_matrix), method = "average")
plot(cave, hang = -1, main = "Average Dendogram", cex = 0.4)

### Cut the dendrogram such that 3 clusters are produced
rect.hclust(cave, k=4, border=2:4)


# Complete
ideal_number_of_clusters(dissimilarity_matrix, "complete", "wss")
ideal_number_of_clusters(dissimilarity_matrix, "complete", "silhouette")
##ideal_number_of_clusters(dissimilarity_matrix, "complete", "gap_stat")

## Perform hierarchical clustering
ccom <- hclust(dist(dissimilarity_matrix), method = "complete")
plot(ccom, hang = -1, main = "Complete Dendogram", cex = 0.4)

### Cut the dendrogram such that 3 clusters are produced
rect.hclust(ccom, k=4, border=2:4)


# Ward.D2
ideal_number_of_clusters(dissimilarity_matrix, "ward.D2", "wss")
ideal_number_of_clusters(dissimilarity_matrix, "ward.D2", "silhouette")
##ideal_number_of_clusters(dissimilarity_matrix, "ward.D2", "gap_stat")

## Perform hierarchical clustering
cwar <- hclust(dist(dissimilarity_matrix), method = "ward.D2")
plot(cwar, hang = -1, main = "Ward Dendogram", cex = 0.4)

### Cut the dendrogram such that 3 clusters are produced
rect.hclust(cwar, k=4, border=2:4)


# Save the dendrograms
dend1 <- as.dendrogram (csin)
dend2 <- as.dendrogram (cave)
dend3 <- as.dendrogram (ccom)
dend4 <- as.dendrogram (cwar)

# Create a correlation matrix
trees <- dendlist("Single"=dend1, "Average"=dend2, "Complete"=dend3, "Ward"=dend4)
baker <- cor.dendlist(trees, method="baker")

# Plot the correlation matrix
corrplot(baker,
         method="circle",
         type="lower",
         tl.col = "black",
         tl.cex = 1.0,
         cl.cex = 1.0,
         addCoef.col = "white",
         number.cex = 0.7,
         col.lim=c(0.0,1),
         col=COL2("RdBu",n = 20)
)
