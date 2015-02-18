\name{correlation.test}
\alias{correlation.test}
\title{Correlation test}
\usage{
  correlation.test(x, y, data, ...)
}
\arguments{
  \item{x}{a numeric variable}

  \item{y}{anoter numeric variable}

  \item{data}{a dataframe}

  \item{...}{passed on to stats::cor.test}
}
\value{
  cor.test object
}
\description{
  Perform a correlation test between two numeric variables
}
\examples{
correlation.test("vs", "mpg", data=mtcars)
}

