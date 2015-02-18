\name{multi.cor}
\alias{multi.cor}
\title{Multiple correlations}
\usage{
  multi.cor(data, vars, type = "pearson")
}
\arguments{
  \item{data}{a data frame}

  \item{vars}{(optional) a vector with the subset of
  variable names}

  \item{type}{type of correlation. Defaults to "pearson".}
}
\value{
  correlation matrix
}
\description{
  Do an analysis with multiple correlations.
}
\examples{
multi.cor(iris)
}

