\name{panel.metaCount}
\alias{panel.metaCount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Panel Function for Counting Events by category
}
\description{
metaCount Panel Function
}

\arguments{
  \item{x}{
Provided by lattice
}
  \item{y}{
Provided by lattice
}
\item{end}{
  Vector containing the ending positions of the events being
  counted. For point events, the vector should contain \code{NA}
  }
  \item{subscripts}{
Provided by lattice
}
  \item{patientid}{
    Data indicating the patient/sample an observation came from
}
\item{scale.factor}{
      Numeric. Indicates maximum height of a bar (for values to be
}
\item{logscale}{
  Logical. If TRUE heights of count bars are in capped log scale. Defaults
    to TRUE
}
\item{logbase}{
      Numeric. Indicates the base of the log to use if logscale =
    TRUE. Ignored if logscale=FALSE. Defaults to 1.5.
}
\item{at.baseline}{
  Logical. If TRUE all count bars for a single category have their base at
    the same point. If FALSE, bars are centered in the area of the plot
    indicating their respective category. Defaults to TRUE
}
\item{colpalette}{
  Color palette to use in the counting portion of the plot. Color is used
    to indicate relative frequency. Should be a character vector
    containing hex color specifications.

}
\item{legend.step}{
  Numeric. Step size used to differentiate between colors in colpalette.
}
  \item{\dots}{
    Unused.
}
}
\details{
This is designed to be a panel function for use inside a lattice
plotting function call.
}
\value{
Not meaningful
}
\author{
Gabriel Becker
}

\keyword{ ~panel }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
