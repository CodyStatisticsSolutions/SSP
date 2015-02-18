\name{data}
\alias{data}
\docType{data}
\title{Test Data
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(data)}
\format{
  A data frame with 150 observations on the following 18 variables.
  \describe{
    \item{\code{ID}}{a numeric vector}
    \item{\code{NormalDependent}}{a numeric vector}
    \item{\code{NotNormalDependent}}{a numeric vector}
    \item{\code{DifferentGroup2}}{a factor with levels \code{Group A} \code{Group B}}
    \item{\code{NotDifferentGroup2}}{a factor with levels \code{Group A} \code{Group B}}
    \item{\code{DifferentGroup3}}{a factor with levels \code{Alpha} \code{Beta} \code{Theta}}
    \item{\code{NotDifferentGroup3}}{a factor with levels \code{Alpha} \code{Beta} \code{Theta}}
    \item{\code{Paired1}}{a numeric vector}
    \item{\code{Paired2Different}}{a numeric vector}
    \item{\code{Paired2NotDifferent}}{a numeric vector}
    \item{\code{Correlate}}{a numeric vector}
    \item{\code{PositiveCorrelate}}{a numeric vector}
    \item{\code{NegativeCorrelate}}{a numeric vector}
    \item{\code{NoCorrelate}}{a numeric vector}
    \item{\code{ChiGroup1}}{a factor with levels \code{Group A} \code{Group B}}
    \item{\code{ChiGroup2Large}}{a factor with levels \code{Alpha} \code{Beta} \code{Theta}}
    \item{\code{ChiGroup2Small}}{a factor with levels \code{Alpha} \code{Beta} \code{Theta}}
    \item{\code{ChiSquareNotSig}}{a factor with levels \code{Bob} \code{George} \code{Hank} \code{Tim}}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(data)
## maybe str(data) ; plot(data) ...
}
\keyword{datasets}
