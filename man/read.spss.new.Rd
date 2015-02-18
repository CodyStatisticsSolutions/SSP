\name{read.spss.new}
\alias{read.spss.new}
\alias{testdata}
\title{Read SPSS (new)}
\usage{
  read.spss.new(file, ..., to.data.frame = FALSE)
}
\arguments{
  \item{file}{a valid spss data file (.sav)}

  \item{...}{arguments passed to read.spss}

  \item{to.data.frame}{(depricated). Will be ignored.}
}
\value{
  data frame
}
\description{
  Special wrapper for read.spss in order to preserve both
  orriginal encodings and factors.
}
\examples{
spssfile <- system.file(package="StatSolutions", "files/1991GS.sav");
mydata <- read.spss.new(spssfile);
}

