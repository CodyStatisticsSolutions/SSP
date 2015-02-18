\name{chi.square}
\alias{chi.square}
\title{Chi Square test}
\usage{
  chi.square(x, y, data)
}
\arguments{
  \item{x}{A string defining a categorial variable that
  exists in the dataset.}

  \item{y}{A string defining another categorial variable
  that exists in the dataset.}

  \item{data}{The dataset}
}
\value{
  A valid JSON string
}
\description{
  Perform a Chi Square test
}
\examples{
chi.square("vs", "am", data=mtcars)
}

