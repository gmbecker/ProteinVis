\name{pfColMap}
\alias{pfColMap}
\docType{data}
\title{
Mapping between PFAM Symbol and Display Color
}
\description{
This data.frame contains an assigned display color (one of 5 distinct
colors) to each distinct PFAM symbol.
}
\usage{data(pfColMap)}
\format{
  A data frame with 10340 observations on the following 2 variables.
  \describe{
    \item{\code{Symbol}}{a factor containing levels which represent the
      symbol for each distinct entry in PFAM
      }
    \item{\code{color}}{a factor with levels for each of the 5 distinct
      display colors
      }
  }
}
\details{
This dataset is used to determine which color any given PFAM domain is
displayed as in the protein structure visualizations.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(pfColMap)
## maybe str(pfColMap) ; plot(pfColMap) ...
}
\keyword{datasets}
