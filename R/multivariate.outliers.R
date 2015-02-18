#' Calculate Mahalanobis distance for outlier removal
#' Required Libraries : plyr
#' 
#' @Calculate Mahalanobis squared distance 
#' @param df - Input dataset as data frame
#'           - For regression, use the columns of the predictors
#' @param whatvar - which variables to be included in the distance calculation    
#'                - default to the entire dataset     
#' @param p - To determine critical chi-square statistic  
#' @param ... - additional arguments to be passed to function mlr 
#' @return output - input datasets with outliers removed (df.outlier.removed.fulldata)
#'                - statement indicating the duplicate column(s) (dupvarlist)

#options(digits=4)

#library(plyr)
#library(jsonlite)

mahdist <- function(df, whatvar=list(), p=0.001, ...) {
  
  # Variables used for distance calculation #
  whatvar.dist <- unique(paste(substitute(whatvar))[-1])
   if (length(whatvar.dist)==0) {
    df.selected <- df 
    whatvar.dist <- names(df)
   } else { 
    df.selected <- df[,whatvar.dist]
   }
  
  # Check for duplicate columns #
  # 1) Indicate which columns are duplicate and to be removed #
  # 2) Continue running the codes with unique columns #
  dupvar <- which(duplicated(as.matrix(df.selected),MARGIN=2))
  if (length(dupvar)>0) {
   dupvarlist <- whatvar.dist[dupvar]
   df.selected <- df.selected[,!names(df.selected) %in% dupvarlist] 
   dupvarlist <- paste("The following duplicate column(s) have been removed to determine the Mahalanobis distance:",
                       paste(dupvarlist,collapse=",")
                       )
   warning(paste(dupvarlist), call.=FALSE)
  } else {dupvarlist <- NULL}

  
  # Incomplete observations to be retained in the dataset #
  # Index on row.names since df.selected and df share the same row names #
  complete <- row.names(df.selected[complete.cases(df.selected),])
  incomplete <- row.names(df.selected)[!row.names(df.selected) %in% complete]

  # Detect outliers among complete observations #
  df.complete <- df.selected[row.names(df.selected) %in% complete,]
  df.center <- colMeans(df.complete,na.rm=TRUE) 
  df.cov <- cov(df.complete,use="complete.obs")
  md <- mahalanobis(x=df.complete, center=df.center, cov=df.cov, ...)  
  
  # Compare to critical chi-square statistic at p level #
  degrees.of.freedom <- ncol(df.complete)
  chisqv <- qchisq(1-p,degrees.of.freedom)
  non.outliers <- row.names(df.complete[which(md<chisqv),])
  retain <- c(non.outliers,incomplete)
  
  # Output datasets # 
  df.outlier.removed.fulldata <- df[row.names(df) %in% retain,]
  df.outlier.removed.fulldata <- arrange(df.outlier.removed.fulldata,as.numeric(row.names(df.outlier.removed.fulldata)))
  return(df.outlier.removed.fulldata)
}

#################################################################################################

# Test #

# Sample data test #
#testdata <- read.csv("G:\\Statistics Solutions\\testdata.csv")
#dim(testdata)
# 150   5

# Testing without missing data #
#data2 <- mahdist(testdata, whatvar=list(NormalDependent, Normal, Outliers1, Paired4))
#toJSON(data2)
#data3 <- mahdist(data2[[1]], whatvar=list(NormalDependent, Normal, Outliers1, Paired4))
#toJSON(data3)

# Testing using a subset of data to be the entire dataset #
#data4 <- testdata[,c("NormalDependent","Normal","Outliers1","Paired4")]
##data5 <- mahdist(data4)
#toJSON(data5)
#data6 <- mahdist(data5[[1]])
#toJSON(data6)

# Testing with missing data #
#data5 <- mahdist(testdata,whatvar=list(NormalDependent, Normal, Outliers1, Paired4,Paired2Different))
#data6 <- mahdist(data5[[1]],whatvar=list(NormalDependent, Normal, Outliers1, Paired4,Paired2Different))
#data7 <- mahdist(data6[[1]],whatvar=list(NormalDependent, Normal, Outliers1, Paired4))
#toJSON(data7)

#data <- read.csv(file=file.choose())
#data2 <- mahdist(data,whatvar=list(NormalDependent, Normal, Outliers1, Paired4,Paired2Different))
#data4 <- mahdist(data2,whatvar=list(NormalDependent, Normal, Outliers1, Paired4,Paired2Different))

