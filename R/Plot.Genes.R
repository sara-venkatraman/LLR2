#' @description
#' Plot the temporal expression trajectories of one or more genes as
#' smooth lines.
#'
#' @param geneData N-by-n numeric dataframe containing the N temporal
#' expression trajectories of length-n to be plotted. Row names should be
#' the names of the genes.
#' @param timePoints Length-n numeric vector of time points corresponding
#' to the times at which the columns of geneData were recorded.
#' @param plotColors Optional vector of colors to use for each line (gene),
#' of the same number of rows as \code{geneData}.
#' @param plotTitle Optional title to place at the top of the plot,
#' centered, or TRUE if the desired title is a list of the gene names.
#' Can include Markdown/HTML tags.
#' @param plotSubtitle An optional title to place under the plot title,
#' centered. Can include Markdown/HTML tags.
#' @param titleSize Numeric font size of the title text. Default 12.
#' @param points Logical indicating whether the observed data points should
#' be plotted atop their corresponding spline interpolants. Default TRUE.
#' @param pointSize Numeric size of observed data points, if \code{points}
#' is TRUE. Default 1.5.
#' @param plotLegend Logical indicating whether a legend should be displayed.
#' Default TRUE.
#' @param legendSize Numeric font size of the legend text. Default 10.
#' @param legendPos Position of legend. Possible values are "top", "bottom",
#' "left", "right". Default "bottom".
#' @param plotGrid Logical indicating whether background grid lines
#' should be drawn. Default TRUE.
#' @param lineLabels Logical indicating whether or not gene names should
#' be printed at the end (right-hand side) of their plotted trajectories.
#' Default FALSE.
#' @param axisLabels Length-2 list with items "x" and "y" specifying x- and
#' y-axis labels. Can include Markdown/HTML tags. Default
#' \code{list(x="Time (hours)", y="Expression")}.
#' @param lineOpacity Numeric opacity of plotted lines. Default 1 (opaque);
#' decrease for transparency.
#'
#' @return A ggplot2 object.
#' @export
Plot.Genes <- function(geneData, timePoints, plotColors, plotTitle="",
                       plotSubtitle="", titleSize=12, points=TRUE, pointSize=1.5,
                       plotLegend=TRUE, legendSize=10, legendPos="bottom",
                       plotGrid=TRUE, lineLabels=FALSE, lineOpacity=1,
                       axisLabels=list(x="Time", y="Expression")) {

  # Set the colors for each line if none are provided
  if(missing(plotColors)) {
    plotColors <- c(brewer.pal(n=8, name="Dark2"), brewer.pal(n=8, name="Set2"), brewer.pal(n=12, name="Paired"))
    plotColors <- rep(plotColors, length.out=nrow(geneData))
  }

  # Adjust opacity of colors
  plotColors <- alpha(plotColors, lineOpacity)

  # Define many time points at which to evaluate each gene's spline interpolant
  interpTimes <- seq(from=timePoints[1], to=tail(timePoints, 1), length.out=500)

  # Evaluate each gene's spline interpolant at each time point. This is the
  # data that will be plotted.
  interpData <- data.frame(interpTimes, matrix(0, nrow=500, ncol=nrow(geneData)))
  interpData[,-1] <- apply(geneData, 1, function(y) splinefun(x=timePoints, y=y, method="natural")(interpTimes))
  colnames(interpData) <- c("time", rownames(geneData))

  # Reshape interpolation data into long format for use with ggplot
  meltInterpData <- melt(interpData, id.var="time")

  # Construct initial plot: smooth lines for each gene
  p <- ggplot(meltInterpData, aes(x=time, y=value, col=variable)) +
    geom_line() + scale_color_manual(values=plotColors) + theme_bw() +
    labs(x=axisLabels$x, y=axisLabels$y) +
    theme(axis.title.x=element_markdown(), axis.title.y=element_markdown())

  # If desired, draw points of specified size at each observed time
  if(points == TRUE) {
    pointData <- data.frame(timePoints, t(geneData))
    colnames(pointData) <- c("time", rownames(geneData))
    meltPointData <- melt(pointData, id.var="time")
    p <- p + geom_point(data=meltPointData, mapping=aes(x=time, y=value, col=variable), size=pointSize)
  }

  # If plot title is desired but not supplied, set it to be a comma-separated
  # list of the gene names
  if(plotTitle == TRUE)
    plotTitle <- paste(rownames(geneData), collapse=", ")

  # Add the plot title, centered and with line height 1.1, with parsing
  p <- p + labs(title=plotTitle) +
    theme(plot.title=element_text(size=titleSize, hjust=0.5)) +
    theme(plot.title=element_markdown(lineheight=1.1))

  # Add the plot subtitle, centered and with line height 1.1, with parsing
  p <- p + labs(subtitle=plotSubtitle) +
    theme(plot.subtitle=element_text(hjust=0.5)) +
    theme(plot.subtitle=element_markdown(lineheight=1.1))

  # Remove grid lines if not desired
  if(plotGrid == FALSE)
    p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

  # If desired, add gene names as labels next to each line
  if(lineLabels == TRUE)
    p <- p + geom_dl(aes(label=variable), method=list(dl.trans(x=x+0.2), "last.points", cex=0.6)) +
    expand_limits(x=tail(timePoints, 1) + 5)

  # If desired, add a legend with specified font size and position
  if(plotLegend == TRUE)
    p <- p + theme(legend.title=element_blank(), legend.position=legendPos,
                   legend.text=element_text(size=legendSize),
                   legend.background=element_rect(size=0.1, linetype="solid", color="black"))
  else
    p <- p + theme(legend.position="none")

  # Return the ggplot object
  p
}
