\name{dataset.convert}
\alias{dataset.convert}
\title{Convert variables in data}
\usage{
  dataset.convert(data, asnumeric = vector(),
    asfactor = vector(), asordered = vector())
}
\arguments{
  \item{data}{dataset}

  \item{asnumeric}{vector of variable names or indices to
  be converted to numbers}

  \item{asfactor}{vector of variable names or indices to be
  converted to factors}

  \item{asordered}{vector of variable names or indices to
  be converted to ordered factors}
}
\value{
  a new dataset
}
\description{
  Special function to convert variables in a dataset from
  nominal to numeric. When used on datasets imported from
  SPSS, it uses orriginal codings.
}

