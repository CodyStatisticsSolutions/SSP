#' Multiple Imputation
#' Required Libraries : mice, VIM, data.table 
#' 
#' @Perform  Multiple Imputation 
#' @param df - data frame that contains all the required variables
#' @param imputed.varlist - incomplete variables to be imputed (each contains missing values)
#' @param pred.varlist - complete or incomplete to be used to impute other variables       
#' @param visseq - order of variables to be imputed - Default "monotone" which means
#'                        starting with the one with the least number of missing values              
#' @param not.for.pred - imputed variables not to be used as predictor but as imputed variables only
#' @param check.mcar - margin plot of checking mcar of the supplied variables
#'                     since it is a xy-plot at least 2 variables are required 
#' @param multilevel - impute missing data under a 2-level linear multilevel model - Default to FALSE
#' @param multilevel.class - Class Variable in 2-level linear multilevel model e.g. school in school/student
#'                           Only one variable is allowed, Cannot have missing values, either factor or numeric
#' @param multilevel.fixed - fixed effects in 2-level linear multilevel model 
#' @param multilevel.random - random effects in 2-level linear multilevel model
#'                            pred.varlist == multilevel.fixed UNION multilevel.random    
#' @param maxit - number of iterations per imputation set - default to 10
#' @param num.of.imputed.sets - number of imputed datasets to be produced - default to 5
#'                              THe final dataset will be a random sample set from all the imputed sets
#' @param mincor - minimum pairwise correlation between two variables to be considered
#'                 a plausible predictor-imputed pair - default to 0.1
#' @param minpuc - The proportion of usable cases measures how many cases with missing data on the target
#'                 variable actually have observed values on the predictor - default to 0.25 (25%)                  
#' @return output - parameters and statistics

#options(digits=4)
library(mice)
#library(VIM)
library(data.table)
#library(jsonlite)

mi.ssp <- function(df, 
                   imputed.varlist = list(),
                   pred.varlist = list(),
                   visseq = list(monotone),
                   not.for.pred = list(),
                   check.mcar = list(),
                   multilevel = FALSE,
                   multilevel.class = list(), 
                   multilevel.fixed = list(), 
                   multilevel.random = list(), 
                   maxit = 10,
                   num.of.imputed.sets = 5,
                   mincor = 0.1,
                   minpuc = 0.25
                   ) {

  # Select appropriate variables to impute #                    
  imputed.varlist <- paste(substitute(imputed.varlist))[-1]
  pred.varlist <- paste(substitute(pred.varlist))[-1]  
  
  # Check missing pattern of the input variables (complete and incomplete) #
  md.check <- md.pattern(df[,c(imputed.varlist, pred.varlist)])
  nr <- as.numeric(attributes(md.check)$dimnames[[1]])
  md.check <- data.frame(md.check)
  names(md.check)[ncol(md.check)] <- "number.of.missing.variables"
  md.check[,"number.of.rows"] <- nr
  md.check <- md.check[-nrow(md.check),]
  md.check[,"percentage"] <- with(md.check,
                                  round(100*number.of.rows/sum(number.of.rows),2))
  
  # Check missing pairs #
  mp.check <- md.pairs(df[,c(imputed.varlist, pred.varlist)])$mm
  mp.check <- data.frame(mp.check)
  mp.check[,"variable"] <- rownames(mp.check)
  rownames(mp.check) <- NULL
  
  # Recommending variables to impute #
  recommend.vars <- data.frame(quickpred(df[,c(imputed.varlist, pred.varlist)], 
                                         mincor = mincor, minpuc = minpuc))
  recommend.vars[,"variables.to.be.imputed"] <- row.names(recommend.vars)
  row.names(recommend.vars) <- NULL
  
  # Graph checking MCAR between variables #
  # If MCAR, the marginal distribution (blue and red boxplots) of the variable will be similar # 
  # i.e. the variable has nothing to do with the missingness of the other variable #  
  check.mcar <- paste(substitute(check.mcar))[-1]
  if (length(check.mcar)>=2)
   {varcomb <- t(combn(check.mcar,2))
     for (p in 1:nrow(varcomb)) {
      marginplot(df[, c(varcomb[p,1], varcomb[p,2])], col = mdc(1:2), cex = 1,
                 cex.lab = 1, cex.numbers = 1, pch = 19, main = "Check MCAR")
     }
   }  
    
  # Imputation order #
  visseq <- paste(substitute(visseq))[-1]
  
  # Variables to be imputed but not to be used a predictors #
  not.for.pred <- paste(substitute(not.for.pred))[-1]
  
  # Multilevel Variables #
  multilevel.class <- paste(substitute(multilevel.class))[-1]
  multilevel.fixed <- paste(substitute(multilevel.fixed))[-1]
  multilevel.random <- paste(substitute(multilevel.random))[-1]
  
  # At least one variable has to have missing values for mice to work #
  if (multilevel == FALSE)
   {
    df.final <- df[,c(imputed.varlist, pred.varlist)]
    # Append the rest of the variables to the final dataset #
    other.varlist <- names(df)[!names(df) %in% c(imputed.varlist, pred.varlist)]
    if (sum(apply(is.na(df.final ), 2, any)) == 0) {
     all.sets <- NULL  
     final.set <- df.final 
    } else {
      ini <- mice(df.final , maxit = 0, pri = FALSE, seed = 100)
      pred <- ini$pred
       if (length(not.for.pred) > 0) {
         pred[,not.for.pred] <- 0 
        } # Some imputed variables will not be used as predictor #
      meth <- ini$meth
      meth[pred.varlist] <- "" # Make sure the incomplete variables used as predictors will not be imputed #
      imp <- mice(df.final, pred = pred, meth = meth, vis = visseq, maxit = maxit, seed = 100, m = num.of.imputed.sets) 
      all.sets <- lapply(1:num.of.imputed.sets, function(x) complete(imp,x))
      all.sets <- lapply(all.sets, function(x) cbind(x, df[, other.varlist]))
      final.set <- cbind(complete(imp), df[, other.varlist])
    }
   } else {
     df.final <- df[,c(multilevel.class, pred.varlist, imputed.varlist)]
     # Append the rest of the variables to the final dataset #
     other.varlist <- names(df)[!names(df) %in% c(imputed.varlist, pred.varlist, multilevel.class)]
     df.final[,"multilevel.class.int"] <- as.integer(df.final[,multilevel.class])
     fixed.random <- vector(mode="numeric", length = length(pred.varlist))
     fixed.random[match(multilevel.fixed, pred.varlist)] <- 1
     fixed.random[match(multilevel.random, pred.varlist)] <- 2
     pred.string <- c(0, fixed.random, rep(0, length(imputed.varlist)), -2)
     ini <- mice(df.final, maxit = 0, pri = FALSE, seed = 100)
     pred <- ini$pred
     # For multilevel model set all predictorMatrix entry to zero first #
     for (vars in pred.varlist) {
       pred[vars,] <- 0
     }
     for (vars in imputed.varlist) {
       pred[vars,] <- pred.string
     }
     meth <- ini$meth
     meth[c(multilevel.class, pred.varlist, c("multilevel.class.int"))] <- "" # Make sure the incomplete variables used as predictors will not be imputed #
     meth[imputed.varlist] <- "2l.norm"
     imp <- mice(df.final, meth = meth, pred = pred,   
                 vis = visseq, maxit = maxit, seed = 100, m = num.of.imputed.sets) 
        
     # All imputed datasets #
     all.sets <- lapply(1:num.of.imputed.sets, 
                        function(x) subset(complete(imp,x), select = -multilevel.class.int))
     all.sets <- lapply(all.sets, function(x) cbind(x, df[, other.varlist]))
     final.set <- cbind(df[, other.varlist], subset(complete(imp), select = -multilevel.class.int))
  }
  
  return(list(md.check, mp.check, recommend.vars, all.sets, final.set))
}

##################################################################

#testdata <- read.csv(file.choose())
#testdata[sample(1:nrow(testdata), 10), "NotDifferentGroup2"] <- NA
#testdata[sample(1:nrow(testdata), 40), "BigGroup"] <- NA

#mi.ssp(testdata, 
#       imputed.varlist = list(Paired2NotDifferent), 
#       pred.varlist = list(NormalDependent, NotNormalDependent, DifferentGroup2, DifferentGroup3, NotDifferentGroup2,
#                           Paired1, Paired2Different, BigGroup)) -> test1

#mi.ssp(testdata, 
#        imputed.varlist = list(Paired.2.Not.Different), 
#        pred.varlist = list(Normal.Dependent, Not.Normal.Dependent, Different.Group.2, Different.Group.3, Not.Different.Group.2,
#                           Paired.1, Paired.2.Different, Big.Group),
#        check.mcar(Paired.2.Not.Different)) -> test1
#test1

#mi.ssp(testdata,
#       imputed.varlist= list(Paired1, Paired2Different, Paired2NotDifferent),
#       pred.varlist = list(NormalDependent, NotNormalDependent, DifferentGroup2, DifferentGroup3, NotDifferentGroup2, Paired1, Paired2Different, BigGroup),
#       check.mcar=list(Paired1, NormalDependent, BigGroup)) -> test2


#mi.ssp(testdata, 
#       imputed.varlist = list(Paired2NotDifferent), 
#       pred.varlist = list(NormalDependent, NotNormalDependent, DifferentGroup2, DifferentGroup3, NotDifferentGroup2,
#                           Paired1, Paired2Different, BigGroup),
#       mincor = 0.25) -> test11



#mi.ssp(testdata, 
#       imputed.varlist = list(Paired2NotDifferent), 
#       pred.varlist = list(NormalDependent, NotNormalDependent, DifferentGroup2, DifferentGroup3, NotDifferentGroup2,
#                           Paired1, Paired2Different, BigGroup),
#       check.mcar = list(Paired2NotDifferent, NormalDependent, BigGroup)) -> test5

#mi.data.frame(test5)
#test5

