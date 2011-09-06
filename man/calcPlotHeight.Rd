\name{calcPlotHeight}
\alias{calcPlotHeight}

\title{
Calculate Appropriate Height for ProteinVis Plot
}
\description{
Convenience function to calculate the total height of a plot based on
the desired height if plotting pfam domains requires only 1 row
}
\usage{
calcPlotHeight(baseheight, type = "metaCount", pfam)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{baseheight}{
numeric. Height of the plot desired if plotting pfam domains requires only one
row (ie there is no overlap)
}
  \item{type}{
character. Indicates whether metaCountStructPlot or proteinStructPlot
will eventually be called.
}
  \item{pfam}{
    data.frame. PFAM data which will be passed to the plotting function
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Height the plot should be taking into account overlapping PFAM domains
}
\author{
Gabriel Becker
}
\seealso{
\code{\link{metaCountStructPlot}} , \code{\link{proteinStructPlot}}
}
\keyword{ ~protein }

