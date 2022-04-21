## Load temporal gene expression data and prior adjacency matrix

geneData <- read.csv("data-raw/GeneData.csv", row.names=1)
priorMatrix <- read.csv("data-raw/PriorMatrix.csv", row.names=1)
colnames(priorMatrix) <- rownames(priorMatrix)

usethis::use_data(geneData, overwrite = TRUE)
usethis::use_data(priorMatrix, overwrite = TRUE)
