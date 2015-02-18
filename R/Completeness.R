#' Identify complete and incomplete variables prior to multiple Imputation
#' It is always important to select variables that causes MAR in the imputation process
#' 
#' @Perform  Identify complete and incomplete variables 
#' @param df - data frame that contains all the variables (complete and incomplete)
#' @param output - list of complete variables and incomplete variables

completeness <- function(df) {
  incomplete.vars <- names(df)[apply(is.na(df), 2, any)]
  complete.vars <- setdiff(names(df), incomplete.vars)
  return(list(complete.vars = complete.vars, incomplete.vars = incomplete.vars))  
}

