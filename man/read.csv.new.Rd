\name{read.csv.new}
\alias{read.csv.new}
\title{Read CSV}
\usage{
  read.csv.new(...)
}
\arguments{
  \item{...}{arguments passed on to read.csv}
}
\value{
  a data frame
}
\description{
  Wrapper for read.csv
}
\examples{
mytemp <- tempfile();
write.csv(iris, mytemp, row.names=FALSE);
newiris <- read.csv.new(mytemp);
stopifnot(identical(iris, newiris));
}

