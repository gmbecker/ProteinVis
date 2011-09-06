\name{metaCountStructPlot}
\alias{metaCountStructPlot}
\alias{proteinStructPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Visualize Protein Structure Information With or Without of Events within Meta-Classes
}
\description{
Visualize protein structure information while optionally plotting
location and frequency of events along the length of the protein.
}
\usage{
metaCountStructPlot(events, catname = "PRIMARY_TISSUE", position =
c("protpos", "protposend"), pfam, pfamLabels = "featureName",  structPred, hydro,
transMem, sigP, xlim, tmposcol = c("start", "end"), main = NULL, simple
= FALSE, at.baseline = TRUE, logscale = TRUE, logbase = 1.5,
scale.factor = 10, colpalette = rev(brewer.pal(11, "RdYlBu")),
legend.step = 0.01, sampleID, key, subtitle, draw = FALSE)

proteinStructPlot(pfam, structPred, hydro, transMem, sigP, xlim,
tmposcol = c("start", "end"), main = NULL, pfamLabels = "featureName",
draw = FALSE)
}

\arguments{
  \item{events}{
    Data.frame containing the data for the events to count (must include
    both position (start and end) and category information)
}
  \item{catname}{
Name of the column in events which specifies category.
}
  \item{position}{
Name of the columns specifying the start and end of the position where
the event occured on the protein.
}
  \item{pfam}{
data.frame containing PFAM domain data, in the form of columns start,
end, as well as label information (see pfamLabels)
}
\item{pfamLabels}{
  Name of the column in pfam which specifies labels for the domains to
  print in the image
  }
  \item{structPred}{

    data.frame containing secondary structure predictions. Should
    contain a column "start" which indicates position on the protein, as
    well as columns "helix", "strand", and "coil" representing estimated
    probabilities/confidence (on a scale of 0 to 10) of the respective
    structures at that point.
}
  \item{hydro}{
    data.frame containing hydrophobicity data. Should contain a column
    "start" indicating position on the protein and "featureValue"
    indicating the estimated hydrophobicity at that position.
}
  \item{transMem}{
    data.frame containing estimated positions of Transmembrane
    Domains. Should contain columns "start", and "end" indicating
    position of the domain, as well as an optional column "featureType"
    indicating the process used to perform the prediction
}
  \item{sigP}{
    data.frame containing the estimated position of the signal
    peptide. Should contain columns "start" and "end" indicating
    position on the protein.
}
  \item{xlim}{
    Numeric. Indicates range of possible amino acid positions. Eg c(1, 857).
}
  \item{tmposcol}{
    Column names in transMem indicating start and end of domain position.
}
  \item{main}{
Main title for the entire visualization.
}
  \item{simple}{
Logical. If TRUE, generates a simpler version of the event count
    plot. Defaults to false.
}
  \item{at.baseline}{
Logical. If TRUE all count bars for a single category have their base at
    the same point. If FALSE, bars are centered in the area of the plot
    indicating their respective category. Defaults to TRUE
}
  \item{logscale}{
Logical. If TRUE heights of count bars are in capped log scale. Defaults
    to TRUE
}
  \item{logbase}{
    Numeric. Indicates the base of the log to use if logscale =
    TRUE. Ignored if logscale=FALSE. Defaults to 1.5.
}
  \item{scale.factor}{
    Numeric. Indicates maximum height of a bar (for values to be
    compared against). Defaults to 10.
}
  \item{colpalette}{
Color palette to use in the counting portion of the plot. Color is used
    to indicate relative frequency. Should be a character vector
    containing hex color specifications.
}
  \item{legend.step}{
    Numeric. Step size used to differentiate between colors in colpalette.
}
  \item{sampleID}{
    Character. Column name in events which indicates unique sample
    identifier. Used to determine total sample size in a given category.
}
  \item{key}{
    Ignored for now.
}
  \item{subtitle}{
    Character. Subtitle for the entire visualization.
  }
  \item{draw}{
    Logical. Should the plot be drawn to the active graphics
    device. Defaults to FALSE. If FALSE a trellis object representing the
    visualization is returned. Currently draw=TRUE does not work.
    }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 If draw=FALSE, a trellis
    object representing the visualization is returned. Currently
    draw=TRUE does not work (see note).
}

\author{
Gabriel Becker
}
\note{
  Currently draw=TRUE does not work. This seems to be due to namespace
  issues in the required latticeExtra.

  To create the plot, simply store the return value of the function
  while draw=FALSE then call print(<result>).
  }
\note{

  Currently, if two estimated TM domains overlap, their union (all
  positions covered by either estimate) is drawn to the image.
  %%  ~~further notes~~
}



\seealso{
  \code{\link{proteinStructPlot}}
  
}

\keyword{ protein }
\keyword{ visualization }% __ONLY ONE__ keyword per line
