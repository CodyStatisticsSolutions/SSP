#' Calculate partial correlation #
#' Required Libraries : corrplot 
#' @Perform Partial correlation
#' @details main function partcor()
#' @param df - data frame containing of all variables, only complete observations are accepted #
#' @param varlist - list of required variables for which partial correlation will be calculated
#' @param method - Type of correlation, Default is Pearson
#'        If "pairwise.complete.obs" is selected, either covariance or correlation matrices may end up not positive semi-definite, 
#'        as well as NA entries if there are no complete pairs for that pair of variables.   
#'        For cov and var, "pairwise.complete.obs" only works with the "pearson" method. 
#' @param use - Type of obs to be used, Default is complete observations 
#' @return output required parameters and statistics

# Loads modified partial correlation function pcor.modified() #
#source("G:\\Statistics Solutions\\Functions\\SS013 - Partial Correlation\\pcor_modified.R")

#options(digits=4)

#library(corrplot)
#library(jsonlite)

# Main Function #

partcor <- function(df, 
                    varlist = list(),
                    method = c("pearson","kendall","spearman"),
                    use = c("complete.obs","pairwise.complete.obs")) {
  
 # Select required variables, if ignored all variables in df will be selected #
 varlist <- paste(substitute(varlist))[-1]
  if (length(varlist)>0) {
   df <- df[,varlist] 
  }
 # Create partial correlation object # 
 pc <- pcor.modified(df, method = method[1], use = use[1]) 
 # Variable Names #
 varnames <- dimnames(pc$estimate)[[1]]
 # Partial Correlations #
 pcorr <- pc$estimate
 print(corrplot(pcorr, method = "ellipse", title = "Partial Correlation Matrix Plot", 
                type = "lower", mar = c(1,0.5,1,0.5), tl.cex = 0.75))
 # p-value - is the correlation significant? 2-sided pvalue derived from Z score#
 pcorr.pvalue <- pc$p.value
 diag(pcorr.pvalue) <- NA
 # Output #
 return(list(vars = varnames, coeffs = pcorr, p = pcorr.pvalue))
}

##################################################################################################################################
################################################## End of Main Function ##########################################################
##################################################################################################################################

# Testing #

# Using test data #
#testdata <- read.csv("G:\\Statistics Solutions\\testdata.csv")
# Continuous variables #
#contvars <- c("NormalDependent","NotNormalDependent","Paired1","Paired2Different","NegativeCorrelate","Outliers1","TypicalLikert1")
#testdata1 <- testdata[,contvars]
#a <- partcor(testdata1) ; fromJSON(toJSON(a))
#a <- partcor(testdata1, varlist=list(NormalDependent,Paired1,Outliers1))
#partcor(testdata1,method="kendall",use="pairwise.complete.obs") # Pairwise is not allowed for Kendall #
#partcor(testdata1,method="spearman",use="pairwise.complete.obs") # Pairwise is not allowed for Spearman #
#partcor(testdata1,method="pearson",use="pairwise.complete.obs") # Pairwise is not allowed for Spearman #
#partcor(testdata1,method="kendall")
#partcor(testdata1,method="spearman", varlist=list(NormalDependent,Paired1,Outliers1))
#toJSON(partcor(testdata1))
#toJSON(partcor(testdata1,method="spearman", varlist=list(NormalDependent,Paired1,Outliers1)))

#fromJSON(toJSON(partcor(data, varlist=list(Normal, NotNormalDependent, Paired1, Outliers1, NegativeCorrelate)), digits=4))
