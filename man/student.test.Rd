\name{student.test}
\alias{student.test}
\title{Student t-test}
\usage{
  student.test(x, y, data, ...)
}
\arguments{
  \item{x}{factor variable with exactly two levels}

  \item{y}{numeric variable}

  \item{data}{data frame}

  \item{...}{arguments pased on to t.test}
}
\value{
  bunch of t-test info
}
\description{
  Compare two groups in a t-test.
}
\examples{
student.test("vs", "mpg", data=mtcars)
}

