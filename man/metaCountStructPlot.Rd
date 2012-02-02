\name{metaCountStructPlot}
\alias{metaCountStructPlot}
\alias{proteinStructPlot}
\alias{createProteinImages}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Visualize Protein Structure Information With or Without of Events within Meta-Classes
}
\description{
Visualize protein structure information while optionally plotting
location and frequency of events along the length of the protein.
}
\usage{
metaCountStructPlot(events, catname = "PRIMARY_TISSUE", requiredCats =
NULL, position =
c("protpos", "protposend"), pfam, pfamLabels = "featureName",  structPred, hydro,
transMem, sigP, xlim, tmposcol = c("start", "end"), main = NULL,
 logscale = TRUE, logbase = 1.506,
scale.factor = 10, colpalette = rev(brewer.pal(11, "RdYlBu")),
legend.step = 0.01, sampleID, subtitle, draw = TRUE, vertGuides = 10,
sequence.counts = NULL)

proteinStructPlot(pfam, structPred, xlim, tmposcol = c("start", "end"),
main = NULL, pfamLabels = "featureName", draw = TRUE, vertGuides = 10)

createProteinImages(events,catname = "PRIMARY_TISSUE",
requiredCats = NULL, position = c("protpos", "protposend"),  pfam,
pfamLabels = "featureName",structPred, hydro, transMem, sigP, xlim,
tmposcol = c("start", "end"), main = NULL, logscale = TRUE,
logbase = 1.506, scale.factor = 10, colpalette =rev(brewer.pal(11,"RdYlBu")),
legend.step = .01, sampleID, subtitle = "Amino Acid Position",
vertGuides = 10, sequence.counts = NULL, metaFileName =
"VariantPlot.svg", structFileName = "StructurePlot.svg", metaHeight,
structHeight, width)

}

\arguments{
  \item{events}{
    Data.frame containing the data for the events to count (must include
    both position (start and end) and category information)
}
  \item{catname}{
Name of the column in events which specifies category.
}
\item{requiredCats}{
  Character vector specifying the categories that must be displayed in
  the graphic (even in the that they contain no observations). They
  should be specified in the order they are to appear starting from the
  bottom of the metaCount section of the plot.. If there are
  categories observed in the data not found in requiredCats, they are
  added to the end of the list in the order they appear. If \code{NULL},
  categories and order are taken from the observed data.
}
  \item{position}{
Name of the columns specifying the start and end of the position where
the event occured on the protein. The column indicated as the end
  position should contain \code{NA} for point events and ending position
  for events spanning more than one position. 
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
  \item{logscale}{
    Logical. If TRUE heights of count bars are in capped log scale. Defaults
    to TRUE
}
  \item{logbase}{
    Numeric. Indicates the base of the log to use if logscale =
    TRUE. Ignored if logscale=FALSE. Defaults to 1.506.
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
  \item{subtitle}{
    Character. Subtitle for the entire visualization. Appears below the image.
  }
  \item{draw}{
    Logical. Should the plot be drawn to the active graphics
    device. Defaults to FALSE. If FALSE a trellis object representing the
    visualization is returned.
  }
  \item{vertGuides}{
    Integer. Number of vertical guidelines to draw in the images. Lines
    will be equally spaced.
  }
  \item{sequence.counts}{
    data.frame or NULL. Data.frame containing the total number of sequenced
    samples. This is necessary because not all samples will contain at
    least one variant, and thus would not appear in
    \code{events}. Should contain columns "category" and
    "count". Categories must match those which appear in \code{events}
    and \code{requiredCats} exactly. Categories with 0 counts but which
    appear in \code{requiredCats} must be included.
  }
  \item{metaFileName}{
    Character. Indicates location to save the svg of the
    metaCountStructPlot image.
  }
  \item{structFileName}{
    Character. Indicates location to save the svg of the
    proteinStructPlot image.
  }
  \item{metaHeight}{
    Numeric. Indicates height desired for variant/structure plot in
    inches.
  }
  \item{structHeight}{
    Numeric. Indicates height desired for structure only plot in inches.
  }
  \item{width}{
    Numeric. Indicates width desired for both plots in inches.
  }
}

\value{
  In the case of proteinStructPlot and metaCountStructPlot a complex
  image is generate based on the inputs. A trellis object is returned
  describing the plot, after being drawn or not as indicated by
  \code{draw}.
  

  In the case of createProteinImages, two interactive svg images are
  created and saved to the harddrive at the locations specified in
  \code{metaFileName} and \code{structFileName}. The return value for
  this function is not meaningful.
  }
\author{
Gabriel Becker
}
\keyword{ protein }
\keyword{ visualization }% __ONLY ONE__ keyword per line
