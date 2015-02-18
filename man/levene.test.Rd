\name{levene.test}
\alias{levene.test}
\title{Levene test}
\usage{
  levene.test(y, data)
}
\arguments{
  \item{y}{a formula}

  \item{data}{a data frame}
}
\value{
  levene test object
}
\description{
  Perform a levene test
}
\examples{
levene.test(Sepal.Width ~ Species, data=iris)
}

