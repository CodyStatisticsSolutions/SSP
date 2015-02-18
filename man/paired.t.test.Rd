\name{paired.t.test}
\alias{paired.t.test}
\title{Paired t-test}
\usage{
  paired.t.test(y1, y2, data, ...)
}
\arguments{
  \item{y1}{variable name 1}

  \item{y2}{variable name 2}

  \item{data}{data frame}

  \item{...}{arguments passed to t.test}
}
\value{
  lots of output
}
\description{
  Do a t-test on paired observations
}
\examples{
paired.t.test("Sepal.Width", "Sepal.Length", iris);
}

