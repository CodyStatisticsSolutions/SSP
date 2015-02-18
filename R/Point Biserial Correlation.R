#' Calculate polyserial (biserial when the factor/categorical variable is dichotomous) and polychoric correlation #
#' Required Libraries : polycor, corrplot 
#' @Perform Polyserial and polychoric correlation
#' @param df - data frame containing the required variables
#' @param varlist - list of required variables for which partial correlation will be calculated 
#' @param ML - Maximum Likelihood is used to estimate correlation. Default is TRUE 
#' @param x - a numeric variable 
#' @param y - an ordered categorical/numeric variable, can be dichotomous
#' @return output required parameters and statistics

#options("scipen"=0, "digits"=4)

#library(polycor)
#library(corrplot)
#library(jsonlite)


# Calculate Polyserial Correlation (plsc) # 
# If y is dichotomous, the result will be point biserial correlation #
# If y is ordered numeric or factor, the result will be polyserial correlation #
plsc <- function(df, x, y) {
 pls.corr <- polyserial(df[,deparse(substitute(x))],df[,deparse(substitute(y))])
 return(pls.corr)
}

# Calculate Polychoric Correlation (plcc) # 
# If both x and y are dichotomous, the result will be tetrachoric correlation #
# If at least one of x and y is ordered numeric or factor, the result will be polychoric correlation #
plcc <- function(df, x, y) {
 pls.corr <- polychor(df[,deparse(substitute(x))],df[,deparse(substitute(y))])
 return(pls.corr)
}

# Heterogeneous Correlation Matrix (allcc) #
# Default to use complete obs and Maximum Likelihood #
# Compute all Pearson, Polyserial and Polychoric correlation when all variables are provided at once #
# The optimization process may fail if two categorical variables are identical #

hetcor.mod <- function(...) {
  do.call("hetcor", list(...)) 
}

allcc <- function(df, varlist = list(), ML=TRUE, ...) {
 # Select required variables, if ignored all variables in df will be selected #
 varlist <- paste(substitute(varlist))[-1]
  if (length(varlist)>0) {
    df <- df[,varlist] 
  }
 # Create correlation matrix of mixed types #
 all.corr <- hetcor(df, ML=ML, ...)
 # Variable Names #
 varnames <- dimnames(all.corr$correlations)[[1]]
 # Correlations #
 corr <- all.corr$correlations
 print(corrplot(corr, method = "ellipse", title = "Correlation Matrix Plot", 
                type = "lower", mar = c(1,0.5,1,0.5), tl.cex = 0.75, order = "original"))
 # Correlation Types (Pearson/Polyserial/Polychoric) #
 corr.type <- all.corr$type
 diag(corr.type) <- "Identity" 
 # Correlaton Standard Error - used to calibrate p value using Z distribution #
 corr.stderr <- all.corr$std.errors
 diag(corr.stderr) <- NA
 # Calibrate 2-sided p values #
 corr.z <- abs(corr/corr.stderr)
 corr.pvalue <- 2*(1-pnorm(corr.z))
 # Output #
 return(list(var = varnames, coeffs = corr, type = corr.type, p = corr.pvalue))
}


###################################################################################################################################
################################################## End of Main Functions ##########################################################
###################################################################################################################################

# Testing #

# Using test data #
#testdata <- read.csv("G:\\Statistics Solutions\\testdata.csv")
#toJSON(plsc(testdata,NormalDependent,ChiGroup1))
#toJSON(plcc(testdata,ChiGroup1,TypicalLikert3))
#a <- allcc(testdata[,1:10])
#toJSON(allcc(testdata[,1:10]))
#fromJSON(toJSON(allcc(testdata[,1:10],varlist=list(NotDifferentGroup2,Paired2Different,Paired1))))

#data <- read.csv(file=file.choose())
#b <- allcc(data2, varlist=list(Normal,NotDifferentGroup2, Paired1, DifferentGroup2, ChiGroup2Large))
#b
#data2 <- nominal.to.ordinal("ChiGroup2Large", order=c("Alpha", "Beta", "Theta"), data)
#citation(package="polycor")

#fromJSON(toJSON(b))
#allcc(data, varlist=list("Normal", "NotDifferentGroup2", "Paired1", "DifferentGroup2", "ChiGroup2Large"))
