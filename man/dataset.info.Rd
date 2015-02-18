\name{dataset.info}
\alias{dataset.info}
\title{Dataset info}
\usage{
  dataset.info(data)
}
\arguments{
  \item{data}{the dataframe}
}
\value{
  list with all kind of summary stuff
}
\description{
  Get some summary and descriptive statistics on a
  dataframe to show in a GUI.
}
\examples{
spssfile <- system.file(package="StatSolutions", "files/1991GS.sav");
mydata <- read.spss.new(spssfile);
dataset.info(mydata);
}

