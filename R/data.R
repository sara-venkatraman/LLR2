#' @title Example temporal gene expression data from Schlamp et al. 2021
#' @description A data set containing gene expression measurements for 1735 genes at 17 time points.
#' @format A data frame with 1735 rows and 17 columns. Row names are the gene names and column names are the time point names.
"geneData"

#' @title Example prior information/adjacency matrix for example dataset
#' @description A symmetric "adjacency" matrix for the example gene expression dataset ('geneData'), as described in Section 3 of Venkatraman et al. Entry (i,j) is 1 if there is prior evidence of association between the i-th and j-th genes, 0 if the association is unlikely, and NA if the association is unknown.
#' @format A symmetric matrix (data frame) with 1735 rows and 1735 columns. Row and column names are the gene names.
"priorMatrix"
