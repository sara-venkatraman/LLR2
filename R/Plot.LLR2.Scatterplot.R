#' @description
#' Draw a scatterplot of lead-lag R^2 (LLR2) values, with LLR2_other on the
#' horizontal axis and LLR2 - LLR2_own on the vertical axis. Examples can
#' be found in Figure 8 of Venkatraman et al. 2021.
#'
#' @param LLR2Mat Length-3 list of N-by-N matrices named \code{LLR2Mat},
#' \code{LLR2Mat.own}, and \code{LLR2Mat.other}. This list can be obtained
#' by running the function \code{LLR2}.
#' @param priorMatrix N-by-N prior adjacency matrix.
#' @param geneSubset Length-k list of gene names, where k \eqn{leq} N,
#' whose pairwise LLR2 values are to be included in the scatterplot.
#' @param interactive Logical indicating whether the plot should be made
#' interactive, i.e. with labels displayed upon hovering.
#' @param plotTitle Optional title to place at the top of the plot,
#' centered. Can include Markdown/HTML tags.
#'
#' @return A ggplot2 object.
#' @export
Plot.LLR2.Scatterplot <- function(LLR2Matrices, priorMatrix, geneSubset,
                                  interactive=TRUE,
                                  plotTitle="Lead-lag R<sup>2</sup> values") {
  # Define a function for vectorizing a symmetric matrix
  Vectorize.Symm.Mat <- function(symmMat) {
    namePairs <- t(combn(colnames(symmMat), 2))
    data.frame(namePairs, value=symmMat[namePairs])
  }

  # Turn LLR2 matrices into pairwise distance lists for both axes
  xAxis <- Vectorize.Symm.Mat(LLR2Matrices$LLR2Mat.other[geneSubset, geneSubset])
  yAxis <- Vectorize.Symm.Mat(LLR2Matrices$LLR2Mat[geneSubset, geneSubset] -
                                LLR2Matrices$LLR2Mat.own[geneSubset, geneSubset])

  # Turn priorMatrix into pairwise prior list, treated as factor
  priorVec <- Vectorize.Symm.Mat(priorMatrix[geneSubset, geneSubset])
  priorVec$value[is.na(priorVec$value)] <- "NA"
  priorVec$value <- factor(priorVec$value, levels=c("0", "NA", "1"))

  # Arrange scatterplot data into a dataframe
  plotData <- data.frame(xAxis=round(xAxis$value, 3),
                         yAxis=round(yAxis$value, 3),
                         prior=priorVec$value)

  # Re-order data so that points with prior=0 are drawn first, and points
  # with prior=1 are drawn last.
  plotData <- plotData[order(plotData$prior),]

  # Create gene name labels for interactive scatterplot
  pointLabels <- paste(xAxis$X1, xAxis$X2, sep=", ")

  # Define gray, blue, and red colors for the scatterplot
  g <- alpha("gray67", 0.9)
  b <- alpha("navy", 0.65)
  r <- alpha("orangered3",0.8)

  # Construct initial plot
  p <- ggplot(plotData, aes(x=xAxis, y=yAxis, color=prior, text=pointLabels)) +
    geom_point(size=0.9) + theme_light() +
    scale_color_manual(name="Prior", labels=c("0", "NA", "1"), values=c(g, b, r)) +
    theme(plot.title=element_text(hjust=0.5)) +
    labs(title=plotTitle, x="LLR<sup>2</sup><sub>other</sub>",
         y="LLR<sup>2</sup> - LLR<sup>2</sup><sub>own</sub>")

  # Further settings for interactive scatterplot
  if(interactive == TRUE) {
    ggplotly(p) %>% layout(legend=list(title=list(text="Prior"), orientation="h",
                                       xanchor="center", x=0.5, y=-0.2))
  }
  # Further settings for non-interactive scatterplot
  else {
    p + theme(plot.title=element_markdown(), axis.title.x = element_markdown(), axis.title.y = element_markdown()) +
      theme(legend.position="bottom", legend.background=element_rect(size=0.1, linetype="solid", color="black")) +
      guides(color=guide_legend(override.aes=list(size=3)))
  }
}
