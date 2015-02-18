\name{linear.regression}
\alias{linear.regression}
\title{Linear Regression}
\usage{
  linear.regression(formula, data, ...)
}
\arguments{
  \item{formula}{regression formula}

  \item{data}{data frame}

  \item{...}{other stuff passed to lm()}
}
\value{
  lots of output
}
\description{
  Perform a linear regression analysis with one or more
  predictors
}
\examples{
linear.regression(speed~dist, data=cars)
}

