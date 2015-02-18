\name{oneway.anova}
\alias{oneway.anova}
\title{Oneway Anova}
\usage{
  oneway.anova(x, y, data, ...)
}
\arguments{
  \item{x}{Variable name 1}

  \item{y}{Variable name 2}

  \item{data}{data frame}

  \item{...}{arguments passed on to lm()}
}
\value{
  lots of output
}
\description{
  Analysis of Variance for a linear model
}
\examples{
oneway.anova("Species", "Sepal.Width", iris)
}

