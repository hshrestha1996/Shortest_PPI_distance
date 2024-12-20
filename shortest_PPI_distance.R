library(data.table)
library(igraph)
getwd()
setwd("./")

PPI <- read.table("./input/Example_PPI_data.sif", header = T)
PPI <- PPI[,c(1,3)]
colnames(PPI) <- c("a", "b")
binaryPPI <- setDT(PPI)


calculatePPIlength <- function(binaryInteractionTable){
  ## get all protein ids in interaction table
  proteins <- unique(c(binaryInteractionTable$a,binaryInteractionTable$b))
  ## create undirected graph using the igraph package
  g <- graph_from_data_frame(binaryInteractionTable, directed = FALSE)
  ## calculate shortest path between all vertices
  distMatrix <- distances(g)
  ## Output as a long format data.table
  combinations <- reshape2::melt(distMatrix)
  combinations <- as.data.table(combinations)
  names(combinations) = c("x","y","dist")
  
  # Ensure dist is numeric and filter out 'Inf' distances
  combinations$dist <- as.numeric(combinations$dist)  # Convert to numeric
  combinations <- combinations[is.finite(dist)]       # Exclude rows where dist is Inf
  
  combinations$x <- as.character(combinations$x)
  combinations$y <- as.character(combinations$y)
  
  #return(combinations)
  # Convert to matrix-like format using dcast
  matrixOutput <- dcast(combinations, x ~ y, value.var = "dist", fill = NA)
  
  return(matrixOutput)
}

pathLength <- calculatePPIlength(binaryPPI)
write.csv(pathLength, "./output/shortest_PPI_distance.csv")
