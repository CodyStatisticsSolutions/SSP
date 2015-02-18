\name{mediation}
\alias{mediation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Moderation
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mediation(iv, mv, dv, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iv}{
%%     ~~Describe \code{iv} here~~
}
  \item{mv}{
%%     ~~Describe \code{mv} here~~
}
  \item{dv}{
%%     ~~Describe \code{dv} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (iv, mv, dv, data) 
{
    output <- list()
    linear1 <- lm(dv ~ iv, data = data, x = TRUE, y = TRUE)
    linear2 <- lm(mv ~ iv, data = data, x = TRUE, y = TRUE)
    linear3 <- lm(dv ~ iv + mv, data = data, x = TRUE, y = TRUE)
    Fstat1 <- as.list(summary(linear1)$fstatistic)
    Fstat2 <- as.list(summary(linear2)$fstatistic)
    Fstat3 <- as.list(summary(linear3)$fstatistic)
    sd.x1 <- apply(linear1$x, 2, sd)
    sd.y1 <- sd(linear1$y)
    std.coef1 <- coef(linear1) * (sd.x1/sd.y1)
    output$ANOVA.table1 <- apply(as.data.frame(anova(linear1)), 
        1, as.list)
    output$COEF.table1 <- apply(cbind(std.coef1, as.data.frame(summary(linear1)$coefficients)), 
        1, as.list)
    sd.x2 <- apply(linear2$x, 2, sd)
    sd.y2 <- sd(linear2$y)
    std.coef2 <- coef(linear2) * (sd.x2/sd.y2)
    output$ANOVA.table2 <- apply(as.data.frame(anova(linear2)), 
        1, as.list)
    output$COEF.table2 <- apply(cbind(std.coef2, as.data.frame(summary(linear2)$coefficients)), 
        1, as.list)
    sd.x3 <- apply(linear3$x, 2, sd)
    sd.y3 <- sd(linear3$y)
    std.coef3 <- coef(linear3) * (sd.x3/sd.y3)
    output$ANOVA.table3 <- apply(as.data.frame(anova(linear3)), 
        1, as.list)
    output$COEF.table3 <- apply(cbind(std.coef3, as.data.frame(summary(linear3)$coefficients)), 
        1, as.list)
    output$model1$p.value <- pf(Fstat1[[1]], Fstat1[[2]], Fstat1[[3]], 
        lower.tail = FALSE)
    output$model2$p.value <- pf(Fstat2[[1]], Fstat2[[2]], Fstat2[[3]], 
        lower.tail = FALSE)
    output$model3$p.value <- pf(Fstat3[[1]], Fstat3[[2]], Fstat3[[3]], 
        lower.tail = FALSE)
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
