\name{read.table.new}
\alias{read.table.new}
\title{Read.table}
\usage{
  read.table.new(...)
}
\arguments{
  \item{...}{arguments passed to read.table}
}
\value{
  data frame
}
\description{
  Wrapper for read.table
}
\examples{
mytemp <- tempfile();
write.table(iris, mytemp, row.names=FALSE);
newiris <- read.table.new(mytemp, header=TRUE);
stopifnot(identical(iris, newiris));
}

