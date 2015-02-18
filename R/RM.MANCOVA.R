#' MANOVA and MANCOVA
#' Required Libraries : car, reshape2, plyr, MVN, biotools, corrplot, lsmeans, data.table 
#' 
#' @Perform  Repeated Measure ANOVA/ANCOVA analysis
#' @param df - data frame that contains all the required variables, in long/vertical format
#' @param .factor1 - Specify levels for the first factor of both 1-way and 2-way RP ANova
#'                    e.g. .factor1 <- list(locationA, locationB, locationC) for a 3-level factor                      
#' @param .factor2 - Specify levels for the second factor of a 2-way RP Anova (can just be test or pretest-posttest, can be higher)
#'                    e.g. .factor2 <- list() for 1-way or .factor2 <- list(pretest,posttest) for 2-way with 2 levels   
#' @param .factor.int - Any interaction involving either .factor1 or .factor2 with indepvars or covars
#'                       mainly for LS means calculation
#' @param depvars - dependent variables, can be repeated measures (within subject variables)
#' @param indepvars - independent variables (between subject variables)
#' @param indepvars.int - interaction terms of indepvars 
#' @param covars - covariates that reduce the error variance
#' @param covars.int - interaction terms of covars / covars:indepvars 
#' @return output - parameters and statistics

#options(digits=4)

#library(car) # Anova #
#library(reshape2) # data transposition #
#library(plyr)
#library(MVN) # Mardia #
#library(biotools) # Box's M #
#library(corrplot)
#library(lattice) # splom #
#library(lsmeans) # marginal means #
#library(data.table)
#library(jsonlite)

##################################################################################################################################
################################################## Start of Main Function ########################################################
##################################################################################################################################

rm.ssp <- function(df, 
                    .factor1 = list(),
                    .factor2 = list(),
                    .factor.int = list(),
                    depvars = list(),
                    indepvars = list(),
                    indepvars.int = list(),
                    covars = list(),
                    covars.int = list()) {
  

  # Convert variable lists into vectors of variable names for assumption test #
  # Make sure every variable specified in the parameters are in df.complete #
  .factor1 <- paste(substitute(.factor1))[-1]
  .factor2 <- paste(substitute(.factor2))[-1]
  .factor.int <- paste(substitute(.factor.int))[-1]
  depvars <- paste(substitute(depvars))[-1]
  indepvars <- paste(substitute(indepvars))[-1]
  indepvars.int <- paste(substitute(indepvars.int))[-1]
  covars <- paste(substitute(covars))[-1]
  covars.int <- paste(substitute(covars.int))[-1]
  
  # 1) Use complete obs #
  # 2) Indicate which columns are duplicate and to be removed #
  # 3) Continue running the codes with unique columns #
  df.complete <- df[complete.cases(df),]
  dupvar <- which(duplicated(as.matrix(df.complete), MARGIN=2))
  if (length(dupvar)>0) {
    df.complete <- df.complete[,!names(df.complete) %in% names(df.complete)[dupvar]]
  }
  
  #####################
  # Check assumptions #
  #####################
  
  # Multivariate outliers will be handled by the Mahalanobis distance module #
  # To remind user that multivariate outliers should be addressed #
  
  # 1) Linearity between DVs (MANOVA/MANCOVA), CVs and DV/CV pairs (MANCOVA) #
  # Identify all continuous variables #
  # Use Scatter Plot Matrices (splom) #
  if (length(covars)==0) {
    contvars <- depvars
    varlabels <- paste0("DV:\n",depvars)
  } else {
    contvars <- c(depvars, covars)
    varlabels <- c(paste0("DV:\n",depvars),paste0("CV:\n",covars))    
  }
  
  print(splom(df.complete[,contvars], 
              main = "Scatter Plot Matrix", 
              xlab = "", axis.text.cex = 0,
              varnames = varlabels))
   
  # 2) Univariate and Multivariate Normality Assumption #
  # Overall instead of by group #
  if (nrow(df.complete)<5000) {
    uvn.final.results <- do.call("rbind.data.frame", lapply(df.complete[,contvars],shapiro.test))[,c("statistic","p.value")]
    names(uvn.final.results)[1] <- "Shapiro-Wilk Statistic"
    uvn.final.results[,"variable"] <- rownames(uvn.final.results)
    uvn.final.results[,"Normally Distributed?"] <-
      with(uvn.final.results,ifelse(p.value<0.05,"NO","YES"))
    rownames(uvn.final.results) <- NULL
  } else {uvn.final.results <- "Sample size is over 5000; The data is assumed to be normally distributed as a result of the Central Limit Theorem"}
  
  # Mardia's Test to examine multivariate normality in each cell #
  # Overall instead of by group #
  mvn.final.results <- mardiaTest(df.complete[,contvars], cov = TRUE, qqplot = FALSE)
  mvn.final.results <- list(skewness.chi.square = mvn.final.results@chi.skew,
                            skewness.p.value = mvn.final.results@p.value.skew,
                            kurtosis.z.statistic = mvn.final.results@z.kurtosis,
                            kurtosis.p.value = mvn.final.results@p.value.kurt)
  
    
  # 3) Homogeneity of Variance of each DV across all IV groups - Bartlett's and Levene's Test #
  # Bartlett's test - can produce boxplots to provide visual inspection #
  if (length(indepvars)>0) {
   bartlett.test.results <- list()
   for (iv in indepvars) {
     bartlett.test.results[[iv]] <- lapply(df.complete[,depvars], bartlett.test, df.complete[,iv])
   }
   bartlett.final.results <- data.frame(matrix(unlist(lapply(bartlett.test.results, function(x) lapply(x, function(y) c(y$statistic, y$parameter, y$p.value)))), ncol = 3, byrow = TRUE), stringsAsFactors = FALSE)  
   names(bartlett.final.results) <- c("Bartlett's K-squared Statistics","df","p.value")
   bartlett.final.results[,"depvar"] <- depvars
   bartlett.final.results[,"indepvar"] <- rep(names(bartlett.test.results), each = length(depvars))
  
   # Levene's test #
   levene.test.results <- list()
   for (iv in indepvars) {
     levene.test.results [[iv]] <- lapply(df.complete[,depvars], leveneTest, as.factor(df.complete[,iv]))
   }
   levene.final.results <- data.frame(matrix(unlist(lapply(levene.test.results, function(x) lapply(x, function(y) c(y$"F value"[1], y$Df[2], y$"Pr(>F)"[1])))), ncol = 3, byrow = TRUE),stringsAsFactors = FALSE)  
   names(levene.final.results) <- c("F-Statistic","df","p.value")
   levene.final.results[,"depvar"] <- depvars
   levene.final.results[,"indepvar"] <- rep(names(levene.test.results), each = length(depvars)) 
  
   # 4) Homogeneity of Covariance of each DV across all IV groups - Box's M-test #
   boxm.test.results <- list()
   for (iv in indepvars) {
     boxm.test.results[[iv]] <- boxM(df.complete[,depvars], df.complete[,iv])
   }
   boxm.final.results <- data.frame(matrix(unlist(lapply(boxm.test.results, function(x) c(x$statistic, x$parameter, x$p.value))),ncol=3,byrow=TRUE),stringsAsFactors=FALSE)  
   names(boxm.final.results) <- c("Chi-square Statistics","df","p.value")
   boxm.final.results[,"depvar"] <- paste(depvars, collapse = ",")
   boxm.final.results[,"indepvar"] <- indepvars
   
   # 7) IV frequency #
   iv.freq <- lapply(as.data.frame(df.complete[,indepvars]), function(x) as.data.frame(table(x)))
   
  } else {
    bartlett.final.results <- NULL
    levene.final.results <- NULL
    boxm.final.results <- NULL
    iv.freq <- NULL
  }
  
  # 5) Homogeneity of regression slopes (between DVs and covariates), if covariates are provided #
  if (length(covars)>0 & length(indepvars)>0) {
   interactions <- apply(expand.grid(indepvars,covars), 1, function(x) paste(x, collapse = "*"))
   rhs <- apply(matrix(interactions,ncol=length(covars)), 1, function(x) paste(x, collapse = "+"))
   models <- apply(expand.grid(depvars,rhs), 1, function(x) paste(x, collapse = "~"))
   models <- paste0("Anova(lm(", models, ", data = df.complete", "), type = 3)")
   models.results <- lapply(models, function(x) eval(parse(text=x)))
   # Add x and y terms to the anova tables #
   terms <- unlist(lapply(models.results, function(x) row.names(x)))
   lhs <- unlist(lapply(models.results, function(x) response <- gsub("Response: ","",attributes(x)$heading[2])))
   lhs <- rep(lhs,unlist(lapply(models.results, function(x) dim(x)[1])))
   # Final anova table #  
   equal.slopes <- ldply(models.results)
   equal.slopes[,"depvar"] <- lhs
   equal.slopes[,"term"] <- terms} else {equal.slopes <- NULL}   
  
  
  # 6) Correlation of DV - Highly correlated DVs weaken the power of the analysis #
  # There should be reasonable level of correlations between DVs #
  dvcorr <- cor(df.complete[,contvars])
  print(corrplot(dvcorr, method = "ellipse", title = "Correlation Matrix Plot", 
                 type = "lower", mar = c(1,0.5,1,0.5), tl.cex = 0.75, order = "original"))
  
  ####################
  # MANOVA / MANCOVA #
  ####################
  
  # Process divided into two parts #
  # 1) Run MANOVA, determine marginal means, multivariate and univariate outcomes #
  # 2) If covariates are provided, run MANCOVA to get the above outcomes #
  
  # Variables to be inserted into the lm formula # 
  depvarlist <- paste(depvars, collapse = ",")
  indepvarlist <- paste(c(indepvars, indepvars.int), collapse = "+")
  
  #######################################
  # Analyses without covariates (ANOVA) #
  #######################################
  
  if (length(.factor2) < 2) 
   {numoffactors <- "one.factor" # One way #
  } else {numoffactors <- "two.factors" } # Two way #
  
  # Set up switch #
  
  phase.switch <- function(phase.type) {
   switch(phase.type,
    one.factor = {
     # Mauchly's Test of Sphericity for within-subject effect #
     factor1 <- factor(.factor1)    
     idata <- data.frame(factor1)       
     null.model <- lm(as.matrix(df.complete[,depvars]) ~ 1)
     sphericity <- Anova(null.model, 
                         idata = idata , 
                         idesign = ~ factor1, 
                         icontrasts = c("contr.sum", "contr.poly"), 
                         type = 3)
     sphericity <- summary(sphericity, multivariate=F)

     # Output sphericity test results #
     Mauchly.test <- data.frame(unclass(sphericity$sphericity.tests))
      names(Mauchly.test)[1] <- "Mauchly's W"
     Mauchly.test[,"Within.Subject.Effect"] <- rownames(Mauchly.test) 
      rownames(Mauchly.test) <- NULL
     GG.HF <- data.frame(unclass(sphericity$pval.adjustments))[,c(1,3)]
     GG.HF[,"Within.Subject.Effect"] <- rownames(GG.HF)
      rownames(GG.HF) <- NULL
     test.of.sphericity <- merge(Mauchly.test, GG.HF, by="Within.Subject.Effect")
     
     # Between and within subjects Effects #
     # Sphericity Assumed and degrees of freedom adjustment #
     if (length(indepvars) == 0 & length(indepvars.int) == 0) {
       manova.results <- eval(parse(text=paste0("lm(cbind(", depvarlist, ") ~ 1, data = df.complete)")))
     } else {
       manova.results <- eval(parse(text=paste0("lm(cbind(", depvarlist, ") ~ ", indepvarlist, ", data = df.complete)"))) 
     }
     aov <- Anova(manova.results, 
                  idata = idata , 
                  idesign = ~ factor1, 
                  icontrasts = c("contr.sum", "contr.poly"), 
                  type = 3)
     aov.summary <- summary(aov, multivariate=F)
     
     # Output effects and p-value adjustment #
     anova.between.within.effects <- data.frame(unclass(aov.summary[[4]]))
     anova.between.within.effects[,"effect"] <- rownames(anova.between.within.effects)
     anova.between.within.effects[,"effect.type"] <- with(anova.between.within.effects,
                                                     ifelse(effect %in% c("(Intercept)", indepvars, indepvars.int), 
                                                      "Between Subject","Within Subject"))
     rownames(anova.between.within.effects) <- NULL
     
     # p-value corrections #
     GG.HF.corrections <- data.frame(aov.summary[[5]])[,c(2,4)]
      names(GG.HF.corrections) <- c("Greenhouse-Geisser Correction","Huynh-Feldt Correction")
     GG.HF.corrections[,"effect"] <- rownames(GG.HF.corrections)
     rownames(GG.HF.corrections) <- NULL
     anova.between.within.effects <- merge(anova.between.within.effects, GG.HF.corrections, by="effect", all.x=TRUE)
          
     # Multivariate Statistics #
     outtests <- car:::print.Anova.mlm
     # allow the function to return the results and disable print
     body(outtests)[[16]] <- quote(invisible(tests))
     body(outtests)[[15]] <- NULL
     # Run the Anova over all tests  
     mtests <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
     manova.mtests.results <- 
      lapply(mtests, function(i) outtests(Anova(manova.results, test.statistic=i, type=3, 
                                                idata = idata , idesign = ~ factor1, 
                                                icontrasts = c("contr.sum", "contr.poly")))[-1,])
     # Stack all results into one data frame #  
     effect <- row.names(manova.mtests.results[[1]])
     multivariate.test <- rep(mtests, lapply(manova.mtests.results, function(x) dim(x)[1]))
     manova.mtests.results <- as.data.frame(do.call("rbind", manova.mtests.results))  
     manova.mtests.results[,"effect"] <- effect
     manova.mtests.results[,"Multivariate test"] <- multivariate.test
     manova.mtests.results[,"effect.type"] <- with(manova.mtests.results,
                                                    ifelse(effect %in% c(indepvars, indepvars.int), 
                                                     "Between Subject","Within Subject"))
     row.names(manova.mtests.results) <- NULL

     # Univariate models to calculate marginal means #
     # Within group #
     # Convert long format data to short format #
     if (length(indepvars.int)>0)
     {all.indepvars <- unlist(strsplit(indepvars.int,split=":"))
      all.indepvars <- unique(c(indepvars, all.indepvars))
     } else {all.indepvars <- indepvars} 
     df.complete.melt <- melt(df.complete[,c(depvars, all.indepvars)], id = c(all.indepvars), measure.vars = depvars, value.name = "depvar")
     df.complete.melt[,"factor1"] <- factor(rep(.factor1, each = nrow(df.complete)), levels = .factor1)

     # LSmeans and posthoc tests #
     if (length(indepvars) == 0 & length(indepvars.int) == 0 & length(.factor.int) == 0) {
       anova.univariate.models <- lm(depvar ~ factor1, data = df.complete.melt)
     } else {
       # Update indepvarlist to incorporate interaction terms involving factors #
       indepvarlist2 <- paste(c(indepvars, indepvars.int, .factor.int), collapse = "+")
       anova.univariate.models <- lm(as.formula(paste("depvar ~ ", indepvarlist2, "+ factor1")), data = df.complete.melt)
     }
     
     # Between and/or Within Subject Marginal means #
     # Default contrast - pairwise #
     # Sample size of each group must be at least 3 #
     anova.univariate.lsmeans.factor1 <- lsmeans(anova.univariate.models, "factor1")
     anova.univariate.lsmeans.factor2 <- anova.univariate.lsmeans.allfactors <- NULL
     anova.univariate.posthoc.factor2 <- NULL
 
     # Post hoc analysis - Pairwise Comparison using Tukey HSD on Marginal means #
     anova.univariate.posthoc.factor1.ci <- as.data.frame(confint(pairs(anova.univariate.lsmeans.factor1)))
     anova.univariate.posthoc.factor1 <- as.data.frame(summary(pairs(anova.univariate.lsmeans.factor1)))
     anova.univariate.posthoc.factor1 <- 
       cbind(anova.univariate.posthoc.factor1,
             anova.univariate.posthoc.factor1.ci[, c("lower.CL","upper.CL")])
     anova.univariate.lsmeans.factor1 <- as.data.frame(summary(anova.univariate.lsmeans.factor1))
     
     # with independent variables #
     if (length(indepvars) > 0) {
     
     # LS means #   
     anova.univariate.lsmeans.allfactors.indepvars <- 
       lapply(indepvars, function(x) lsmeans(anova.univariate.models,c(x,"factor1")))
     
     # Posthoc #
     anova.univariate.lsmeans.allfactors.indepvars <- lapply(anova.univariate.lsmeans.allfactors.indepvars, function(x) summary(x))
     grpname <- rep(indepvars,lapply(anova.univariate.lsmeans.allfactors.indepvars, function(x) dim(x)[1]))
     anova.univariate.lsmeans.allfactors.indepvars <- as.data.frame(rbindlist(anova.univariate.lsmeans.allfactors.indepvars))
     anova.univariate.lsmeans.allfactors.indepvars[,"indepvar"] <- grpname
     names(anova.univariate.lsmeans.allfactors.indepvars)[1] <- "indepvarvalue"
     
     } else {anova.univariate.lsmeans.allfactors.indepvars <- NULL}
     
     
     ###############################
     # RP with covariates (ANCOVA) #
     ###############################

     if (length(covars)==0 & length(covars.int)==0) {
       ancova.between.within.effects <- mancova.mtests.results <-
         ancova.univariate.lsmeans.factor1 <- ancova.univariate.lsmeans.factor2 <- 
         ancova.univariate.posthoc.factor1 <- ancova.univariate.posthoc.factor2 <- 
         ancova.univariate.lsmeans.allfactors <- ancova.univariate.lsmeans.allfactors.indepvars <- NULL
     } else { 
       covarlist <- paste(c(covars, covars.int), collapse="+")
        
       # Between and within subjects Effects #
       # Sphericity Assumed and degrees of freedom adjustment #
       if (length(indepvars) == 0 & length(indepvars.int) == 0) {
         mancova.results <- eval(parse(text = paste0("lm(cbind(", depvarlist, ") ~ ", covarlist, ", data = df.complete)")))
       } else {
         mancova.results <- eval(parse(text = paste0("lm(cbind(", depvarlist, ") ~ ", indepvarlist, " + ", covarlist, ", data = df.complete)")))
       }
       aov <- Anova(mancova.results, 
                    idata = idata , 
                    idesign = ~ factor1, 
                    icontrasts = c("contr.sum", "contr.poly"), 
                    type = 3)
       aov.summary <- summary(aov, multivariate=F)
       
       # Output effects and p-value adjustment #
       ancova.between.within.effects <- data.frame(unclass(aov.summary[[4]]))
       ancova.between.within.effects[,"effect"] <- rownames(ancova.between.within.effects)
       ancova.between.within.effects[,"effect.type"] <- with(ancova.between.within.effects,
                                                            ifelse(effect %in% c("(Intercept)", indepvars, indepvars.int, covars, covars.int), 
                                                                   "Between Subject","Within Subject"))
       rownames(ancova.between.within.effects) <- NULL
       
       # p-value corrections #
       GG.HF.corrections <- data.frame(aov.summary[[5]])[,c(2,4)]
       names(GG.HF.corrections) <- c("Greenhouse-Geisser Correction","Huynh-Feldt Correction")
       GG.HF.corrections[,"effect"] <- rownames(GG.HF.corrections)
       rownames(GG.HF.corrections) <- NULL
       ancova.between.within.effects <- merge(ancova.between.within.effects, GG.HF.corrections, by="effect", all.x=TRUE)
       
       # Multivariate Statistics #
       outtests <- car:::print.Anova.mlm
       # allow the function to return the results and disable print
       body(outtests)[[16]] <- quote(invisible(tests))
       body(outtests)[[15]] <- NULL
       # Run the Anova over all tests  
       mtests <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
       mancova.mtests.results <- 
         lapply(mtests, function(i) outtests(Anova(mancova.results, test.statistic=i, type=3, 
                                                   idata = idata , idesign = ~ factor1, 
                                                   icontrasts = c("contr.sum", "contr.poly")))[-1,])
       # Stack all results into one data frame #  
       effect <- row.names(mancova.mtests.results[[1]])
       multivariate.test <- rep(mtests, lapply(mancova.mtests.results, function(x) dim(x)[1]))
       mancova.mtests.results <- as.data.frame(do.call("rbind", mancova.mtests.results))  
       mancova.mtests.results[,"effect"] <- effect
       mancova.mtests.results[,"Multivariate test"] <- multivariate.test
       mancova.mtests.results[,"effect.type"] <- with(mancova.mtests.results,
                                                     ifelse(effect %in% c(indepvars, indepvars.int, covars, covars.int), 
                                                            "Between Subject","Within Subject"))
       row.names(mancova.mtests.results) <- NULL
       
       # Univariate models to calculate marginal means #
       # Within group #
       # Convert long format data to short format #
       if (length(indepvars.int)>0)
       {all.indepvars <- unlist(strsplit(indepvars.int,split=":"))
        all.indepvars <- unique(c(indepvars, all.indepvars))
       } else {all.indepvars <- indepvars} 
       if (length(covars.int)>0)
       {all.covars <- unlist(strsplit(covars.int,split=":"))
        all.covars <- unique(c(covars, all.covars))
       } else {all.covars <- covars} 
       
       df.complete.melt <- melt(df.complete[,c(depvars, all.indepvars, all.covars)], id = c(all.indepvars, all.covars), measure.vars = depvars, value.name = "depvar")
       df.complete.melt[,"factor1"] <- factor(rep(.factor1, each = nrow(df.complete)), levels = .factor1)
       
       # LSmeans and posthoc tests #
       if (length(indepvars) == 0 & length(indepvars.int) == 0 & length(.factor.int) == 0) {
         ancova.univariate.models <- lm(as.formula(paste("depvar ~ ", covarlist, "+ factor1")), data = df.complete.melt)
       } else {
         # Update indepvarlist to incorporate interaction terms involving factors #
         indepvarlist2 <- paste(c(indepvars, indepvars.int, .factor.int), collapse = "+")
         ancova.univariate.models <- lm(as.formula(paste("depvar ~ ", indepvarlist2, "+", covarlist, "+ factor1")), data = df.complete.melt)
       }
       
       # Between and/or Within Subject Marginal means #
       # Default contrast - pairwise #
       # Sample size of each group must be at least 3 #
       ancova.univariate.lsmeans.factor1 <- lsmeans(ancova.univariate.models, "factor1")
       ancova.univariate.lsmeans.factor2 <- ancova.univariate.lsmeans.allfactors <- NULL
       ancova.univariate.posthoc.factor2 <- NULL
       
       # Post hoc analysis - Pairwise Comparison using Tukey HSD on Marginal means #
       ancova.univariate.posthoc.factor1.ci <- as.data.frame(confint(pairs(ancova.univariate.lsmeans.factor1)))
       ancova.univariate.posthoc.factor1 <- as.data.frame(summary(pairs(ancova.univariate.lsmeans.factor1)))
       ancova.univariate.posthoc.factor1 <- 
         cbind(ancova.univariate.posthoc.factor1,
               ancova.univariate.posthoc.factor1.ci[, c("lower.CL","upper.CL")])
       ancova.univariate.lsmeans.factor1 <- as.data.frame(summary(ancova.univariate.lsmeans.factor1))
       
       # with independent variables #
       if (length(indepvars) > 0) {
         
         # LS means #   
         ancova.univariate.lsmeans.allfactors.indepvars <- 
           lapply(indepvars, function(x) lsmeans(ancova.univariate.models,c(x,"factor1")))
         
         # Posthoc #
         ancova.univariate.lsmeans.allfactors.indepvars <- lapply(ancova.univariate.lsmeans.allfactors.indepvars, function(x) summary(x))
         grpname <- rep(indepvars,lapply(ancova.univariate.lsmeans.allfactors.indepvars, function(x) dim(x)[1]))
         ancova.univariate.lsmeans.allfactors.indepvars <- as.data.frame(rbindlist(ancova.univariate.lsmeans.allfactors.indepvars))
         ancova.univariate.lsmeans.allfactors.indepvars[,"indepvar"] <- grpname
         names(ancova.univariate.lsmeans.allfactors.indepvars)[1] <- "indepvarvalue"
         
       } else {ancova.univariate.lsmeans.allfactors.indepvars <- NULL}     
     }
    },
    
    two.factors = {
     # Mauchly's Test of Sphericity for within subject effect #
     factor2 <- factor(rep(.factor2, each = length(.factor1)), levels = .factor2)
     factor1 <- factor(.factor1)  
     idata <- data.frame(factor2, factor1)
     null.model <- lm(as.matrix(df.complete[,depvars]) ~ 1)
     sphericity <- Anova(null.model, 
                         idata = idata , 
                         idesign = ~ factor2 * factor1, 
                         icontrasts = c("contr.sum", "contr.poly"), 
                         type = 3) 
     sphericity <- summary(sphericity, multivariate=F)
     
     # Output sphericity test results #
     Mauchly.test <- data.frame(unclass(sphericity$sphericity.tests))
     names(Mauchly.test)[1] <- "Mauchly's W"
     Mauchly.test[,"Within.Subject.Effect"] <- rownames(Mauchly.test) 
     rownames(Mauchly.test) <- NULL
     GG.HF <- data.frame(unclass(sphericity$pval.adjustments))[,c(1,3)]
     GG.HF[,"Within.Subject.Effect"] <- rownames(GG.HF)
     rownames(GG.HF) <- NULL
     test.of.sphericity <- merge(Mauchly.test, GG.HF, by="Within.Subject.Effect")
      
     # Between and within subjects Effects #     
     # Sphericity Assumed and degrees of freedom adjustment #  
     if (length(indepvars) == 0 & length(indepvars.int) == 0) {
       manova.results <- eval(parse(text=paste0("lm(cbind(", depvarlist, ") ~ 1, data = df.complete)")))
     } else {
       manova.results <- eval(parse(text=paste0("lm(cbind(", depvarlist, ") ~ ", indepvarlist, ", data = df.complete)"))) 
     }
     aov <- Anova(manova.results, 
                  idata = idata , 
                  idesign = ~ factor2 * factor1, 
                  icontrasts = c("contr.sum", "contr.poly"), 
                  type = 3) 
     aov.summary <- summary(aov)
     
     # Output effects and p-value adjustment #
     anova.between.within.effects <- data.frame(unclass(aov.summary[[4]]))
     anova.between.within.effects[,"effect"] <- rownames(anova.between.within.effects)
     anova.between.within.effects[,"effect.type"] <- with(anova.between.within.effects,
                                                     ifelse(effect %in% c("(Intercept)", indepvars, indepvars.int), 
                                                      "Between Subject","Within Subject"))
     rownames(anova.between.within.effects) <- NULL
     
     # p-value corrections #
     GG.HF.corrections <- data.frame(aov.summary[[5]])[,c(2,4)]
     names(GG.HF.corrections) <- c("Greenhouse-Geisser Correction","Huynh-Feldt Correction")
     GG.HF.corrections[,"effect"] <- rownames(GG.HF.corrections)
     rownames(GG.HF.corrections) <- NULL
     anova.between.within.effects <- merge(anova.between.within.effects, GG.HF.corrections, by="effect", all.x=TRUE)
     
     # Multivariate Statistics #
     outtests <- car:::print.Anova.mlm
     # allow the function to return the results and disable print
     body(outtests)[[16]] <- quote(invisible(tests))
     body(outtests)[[15]] <- NULL
     # Run the Anova over all tests  
     mtests <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
     manova.mtests.results <- 
       lapply(mtests, function(i) outtests(Anova(manova.results, test.statistic=i, type=3, 
                                                 idata = idata , idesign = ~ factor2 * factor1, 
                                                 icontrasts = c("contr.sum", "contr.poly")))[-1,])
     # Stack all results into one data frame #
     effect <- row.names(manova.mtests.results[[1]])
     multivariate.test <- rep(mtests, lapply(manova.mtests.results, function(x) dim(x)[1]))
     manova.mtests.results <- as.data.frame(do.call("rbind", manova.mtests.results))  
     manova.mtests.results[,"effect"] <- effect
     manova.mtests.results[,"Multivariate test"] <- multivariate.test
     manova.mtests.results[,"effect.type"] <- with(manova.mtests.results,
                                                   ifelse(effect %in% c(indepvars, indepvars.int), 
                                                          "Between Subject","Within Subject"))
     row.names(manova.mtests.results) <- NULL
     
     # Univariate models to calculate marginal means #
     # Within group #
     # Convert long format data to short format #
     if (length(indepvars.int)>0)
      {all.indepvars <- unlist(strsplit(indepvars.int,split=":"))
       all.indepvars <- unique(c(indepvars, all.indepvars))
     } else {all.indepvars <- indepvars}  
     df.complete.melt <- melt(df.complete[,c(depvars, all.indepvars)], id = c(all.indepvars), measure.vars = depvars, value.name = "depvar")
     df.complete.melt[,"factor2"] <- factor(rep(.factor2, each = nrow(df.complete) * length(.factor1)), levels = .factor2)
     df.complete.melt[,"factor1"] <- factor(rep(rep(.factor1, each = nrow(df.complete)), length(.factor2)), level = .factor1)

     # LSmeans and posthoc tests #
     if (length(indepvars) == 0 & length(indepvars.int) == 0 & length(.factor.int) == 0) {
       anova.univariate.models <- lm(depvar ~ factor1 * factor2, data = df.complete.melt)
     } else {
       # Update indepvarlist to incorporate interaction terms involving factors #
       indepvarlist2 <- paste(c(indepvars, indepvars.int, .factor.int), collapse = "+")
       anova.univariate.models <- lm(as.formula(paste("depvar ~ ", indepvarlist2, "+ factor1 * factor2")), data = df.complete.melt)
     }     
     # Between and/or Within Subject Marginal means #
     # Default contrast - pairwise #
     # Sample size of each group must be at least 3 #
     anova.univariate.lsmeans.factor2 <- lsmeans(anova.univariate.models,c("factor2"))
     anova.univariate.lsmeans.factor1 <- lsmeans(anova.univariate.models,c("factor1"))
     anova.univariate.lsmeans.allfactors <- lsmeans(anova.univariate.models,c("factor2","factor1"))
     
     # Post hoc analysis - Pairwise Comparison using Tukey HSD on Marginal means #
     # Only report between, timepoints and phase #
     
     # Between factor2 #
     anova.univariate.posthoc.factor2.ci <- as.data.frame(confint(pairs(anova.univariate.lsmeans.factor2)))
     anova.univariate.posthoc.factor2 <- as.data.frame(summary(pairs(anova.univariate.lsmeans.factor2)))
     anova.univariate.posthoc.factor2 <- 
       cbind(anova.univariate.posthoc.factor2,
             anova.univariate.posthoc.factor2.ci[, c("lower.CL","upper.CL")])
     anova.univariate.lsmeans.factor2 <- as.data.frame(summary(anova.univariate.lsmeans.factor2))
     
     # Between factor1 #
     anova.univariate.posthoc.factor1.ci <- as.data.frame(confint(pairs(anova.univariate.lsmeans.factor1)))
     anova.univariate.posthoc.factor1 <- as.data.frame(summary(pairs(anova.univariate.lsmeans.factor1)))
     anova.univariate.posthoc.factor1 <- 
       cbind(anova.univariate.posthoc.factor1,
             anova.univariate.posthoc.factor1.ci[, c("lower.CL","upper.CL")])
     anova.univariate.lsmeans.factor1 <- as.data.frame(summary(anova.univariate.lsmeans.factor1))
     
     # Convert all.factors and between.within to data.frame #
     anova.univariate.lsmeans.allfactors <- as.data.frame(summary(anova.univariate.lsmeans.allfactors))
     
     if (length(indepvars) > 0) {
       
       # LS means #   
       anova.univariate.lsmeans.allfactors.indepvars <- 
         lapply(indepvars, function(x) lsmeans(anova.univariate.models,c(x,"factor2","factor1")))
       
       # Posthoc #
       anova.univariate.lsmeans.allfactors.indepvars <- lapply(anova.univariate.lsmeans.allfactors.indepvars, function(x) summary(x))
       grpname <- rep(indepvars,lapply(anova.univariate.lsmeans.allfactors.indepvars, function(x) dim(x)[1]))
       anova.univariate.lsmeans.allfactors.indepvars <- as.data.frame(rbindlist(anova.univariate.lsmeans.allfactors.indepvars))
       anova.univariate.lsmeans.allfactors.indepvars[,"indepvar"] <- grpname
       names(anova.univariate.lsmeans.allfactors.indepvars)[1] <- "indepvarvalue"
       
     } else {anova.univariate.lsmeans.allfactors.indepvars <- NULL}
     
     ###############################
     # RP with covariates (ANCOVA) #
     ###############################
     
     if (length(covars)==0 & length(covars.int)==0) {
       ancova.between.within.effects <- mancova.mtests.results <-
         ancova.univariate.lsmeans.factor1 <- ancova.univariate.lsmeans.factor2 <- 
         ancova.univariate.posthoc.factor1 <- ancova.univariate.posthoc.factor2 <- 
         ancova.univariate.lsmeans.allfactors <- ancova.univariate.lsmeans.allfactors.indepvars <- NULL
     } else { 
       covarlist <- paste(c(covars, covars.int), collapse="+")
       
       # Between and within subjects Effects #
       # Sphericity Assumed and degrees of freedom adjustment #
       if (length(indepvars) == 0 & length(indepvars.int) == 0) {
         mancova.results <- eval(parse(text = paste0("lm(cbind(", depvarlist, ") ~ ", covarlist, ", data = df.complete)")))
       } else {
         mancova.results <- eval(parse(text = paste0("lm(cbind(", depvarlist, ") ~ ", indepvarlist, " + ", covarlist, ", data = df.complete)")))
       }
       aov <- Anova(mancova.results, 
                    idata = idata , 
                    idesign = ~ factor2 * factor1, 
                    icontrasts = c("contr.sum", "contr.poly"), 
                    type = 3)
       aov.summary <- summary(aov, multivariate=F)
       
       # Output effects and p-value adjustment #
       ancova.between.within.effects <- data.frame(unclass(aov.summary[[4]]))
       ancova.between.within.effects[,"effect"] <- rownames(ancova.between.within.effects)
       ancova.between.within.effects[,"effect.type"] <- with(ancova.between.within.effects,
                                                             ifelse(effect %in% c("(Intercept)", indepvars, indepvars.int, covars, covars.int), 
                                                                    "Between Subject","Within Subject"))
       rownames(ancova.between.within.effects) <- NULL
       
       # p-value corrections #
       GG.HF.corrections <- data.frame(aov.summary[[5]])[,c(2,4)]
       names(GG.HF.corrections) <- c("Greenhouse-Geisser Correction","Huynh-Feldt Correction")
       GG.HF.corrections[,"effect"] <- rownames(GG.HF.corrections)
       rownames(GG.HF.corrections) <- NULL
       ancova.between.within.effects <- merge(ancova.between.within.effects, GG.HF.corrections, by="effect", all.x=TRUE)
       
       # Multivariate Statistics #
       outtests <- car:::print.Anova.mlm
       # allow the function to return the results and disable print
       body(outtests)[[16]] <- quote(invisible(tests))
       body(outtests)[[15]] <- NULL
       # Run the Anova over all tests  
       mtests <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
       mancova.mtests.results <- 
         lapply(mtests, function(i) outtests(Anova(mancova.results, test.statistic=i, type=3, 
                                                   idata = idata , idesign = ~ factor2 * factor1, 
                                                   icontrasts = c("contr.sum", "contr.poly")))[-1,])
       # Stack all results into one data frame #  
       effect <- row.names(mancova.mtests.results[[1]])
       multivariate.test <- rep(mtests, lapply(mancova.mtests.results, function(x) dim(x)[1]))
       mancova.mtests.results <- as.data.frame(do.call("rbind", mancova.mtests.results))  
       mancova.mtests.results[,"effect"] <- effect
       mancova.mtests.results[,"Multivariate test"] <- multivariate.test
       mancova.mtests.results[,"effect.type"] <- with(mancova.mtests.results,
                                                      ifelse(effect %in% c(indepvars, indepvars.int, covars, covars.int), 
                                                             "Between Subject","Within Subject"))
       row.names(mancova.mtests.results) <- NULL
       
       # Univariate models to calculate marginal means #
       # Within group #
       # Convert long format data to short format #
       if (length(indepvars.int)>0)
       {all.indepvars <- unlist(strsplit(indepvars.int,split=":"))
        all.indepvars <- unique(c(indepvars, all.indepvars))
       } else {all.indepvars <- indepvars} 
       if (length(covars.int)>0)
       {all.covars <- unlist(strsplit(covars.int,split=":"))
        all.covars <- unique(c(covars, all.covars))
       } else {all.covars <- covars} 
       
       df.complete.melt <- melt(df.complete[,c(depvars, all.indepvars, all.covars)], id = c(all.indepvars, all.covars), measure.vars = depvars, value.name = "depvar")
       df.complete.melt[,"factor2"] <- factor(rep(.factor2, each = nrow(df.complete) * length(.factor1)), levels = .factor2)
       df.complete.melt[,"factor1"] <- factor(rep(rep(.factor1, each = nrow(df.complete)), length(.factor2)), level = .factor1)
       
       # LSmeans and posthoc tests #
       if (length(indepvars) == 0 & length(indepvars.int) == 0 & length(.factor.int) == 0) {
         ancova.univariate.models <- lm(as.formula(paste("depvar ~ ", covarlist, "+ factor1 * factor2")), data = df.complete.melt)
       } else {
         indepvarlist2 <- paste(c(indepvars, indepvars.int, .factor.int), collapse = "+")
         ancova.univariate.models <- lm(as.formula(paste("depvar ~ ", indepvarlist2, "+", covarlist, "+ factor1 * factor2")), data = df.complete.melt)
       }
      
       # Between and/or Within Subject Marginal means #
       # Default contrast - pairwise #
       # Sample size of each group must be at least 3 #
       ancova.univariate.lsmeans.factor2 <- lsmeans(ancova.univariate.models,c("factor2"))
       ancova.univariate.lsmeans.factor1 <- lsmeans(ancova.univariate.models,c("factor1"))
       ancova.univariate.lsmeans.allfactors <- lsmeans(ancova.univariate.models,c("factor2","factor1"))
       
       # Post hoc analysis - Pairwise Comparison using Tukey HSD on Marginal means #
       # Only report between, timepoints and phase #
       
       # Between factor2 #
       ancova.univariate.posthoc.factor2.ci <- as.data.frame(confint(pairs(ancova.univariate.lsmeans.factor2)))
       ancova.univariate.posthoc.factor2 <- as.data.frame(summary(pairs(ancova.univariate.lsmeans.factor2)))
       ancova.univariate.posthoc.factor2 <- 
         cbind(ancova.univariate.posthoc.factor2,
               ancova.univariate.posthoc.factor2.ci[, c("lower.CL","upper.CL")])
       ancova.univariate.lsmeans.factor2 <- as.data.frame(summary(ancova.univariate.lsmeans.factor2))
       
       # Between factor1 #
       ancova.univariate.posthoc.factor1.ci <- as.data.frame(confint(pairs(ancova.univariate.lsmeans.factor1)))
       ancova.univariate.posthoc.factor1 <- as.data.frame(summary(pairs(ancova.univariate.lsmeans.factor1)))
       ancova.univariate.posthoc.factor1 <- 
         cbind(ancova.univariate.posthoc.factor1,
               ancova.univariate.posthoc.factor1.ci[, c("lower.CL","upper.CL")])
       ancova.univariate.lsmeans.factor1 <- as.data.frame(summary(ancova.univariate.lsmeans.factor1))
       
       # Convert all.factors and between.within to data.frame #
       ancova.univariate.lsmeans.allfactors <- as.data.frame(summary(ancova.univariate.lsmeans.allfactors))
       
       if (length(indepvars) > 0) {
         
         # LS means #   
         ancova.univariate.lsmeans.allfactors.indepvars <- 
           lapply(indepvars, function(x) lsmeans(ancova.univariate.models,c(x,"factor2","factor1")))
         
         # Posthoc #
         ancova.univariate.lsmeans.allfactors.indepvars <- lapply(ancova.univariate.lsmeans.allfactors.indepvars, function(x) summary(x))
         grpname <- rep(indepvars,lapply(ancova.univariate.lsmeans.allfactors.indepvars, function(x) dim(x)[1]))
         ancova.univariate.lsmeans.allfactors.indepvars <- as.data.frame(rbindlist(ancova.univariate.lsmeans.allfactors.indepvars))
         ancova.univariate.lsmeans.allfactors.indepvars[,"indepvar"] <- grpname
         names(ancova.univariate.lsmeans.allfactors.indepvars)[1] <- "indepvarvalue"
         
       } else {ancova.univariate.lsmeans.allfactors.indepvars <- NULL}
    
     }
     
    }
   )
   
   return(list(test.of.sphericity, 
               anova.between.within.effects, manova.mtests.results,
               anova.univariate.lsmeans.factor1, anova.univariate.lsmeans.factor2,
               anova.univariate.lsmeans.allfactors, anova.univariate.lsmeans.allfactors.indepvars,
               anova.univariate.posthoc.factor1, anova.univariate.posthoc.factor2,
               ancova.between.within.effects, mancova.mtests.results,
               ancova.univariate.lsmeans.factor1, ancova.univariate.lsmeans.factor2,
               ancova.univariate.lsmeans.allfactors, ancova.univariate.lsmeans.allfactors.indepvars,
               ancova.univariate.posthoc.factor1, ancova.univariate.posthoc.factor2
               ))
  }  
  
  # Output Mauchly's W, between and within subject effects and multivariate statistics #
  all.effects <- phase.switch(numoffactors)  
  
  # Final Output #  
  return(list(uvn = uvn.final.results, mvn = mvn.final.results, 
              bartlett = bartlett.final.results, levene = levene.final.results, 
              boxm = boxm.final.results, slopes = equal.slopes, cor = dvcorr, n = iv.freq, effects = all.effects)) 

}
  

##################################################################################################################################
################################################## End of Main Function ##########################################################
##################################################################################################################################

# Testing #
# Create a test data set in long format #
# 200 subject ids#
# dv1 & dv2 - 5 obs per subject with a total of 1000 rows #
# Take out 100 rows as some subjects do not have complete observations #

#set.seed(100)
#testdata <- data.frame(dv1=rnorm(1000, 100, 30),dv2=rnorm(1000, 100, 20), dv3=rnorm(1000, 100, 20))
#testdata[,"id"] <- rep(1:200, each=5)
#testdata[,"time"] <- rep(1:5, 200)
#testdata <- melt(testdata, id = c("id", "time"))
#testdata <- testdata[-sample(nrow(testdata), size = 50, replace = FALSE),]
#testdata <- dcast(testdata, id ~ variable + time, value.var = "value")

# Add Independent variables iv1 & iv2, and covariates cv1, cv2 & cv3 #
# Some IVs and CVs may have missing values too 

#othvars <- data.frame(id=1:200)
#othvars[,"iv1"] <- as.factor(sample(c("x1","x2","x3"), size = 200, replace = TRUE))
#othvars[,"iv2"] <- as.factor(sample(c("y1","y2","y3","y4"), size = 200, replace = TRUE))
#othvars[,"cv1"] <- rnorm(200, 50, 10)
#othvars[,"cv2"] <- rnorm(200, 60, 15)
#othvars[,"cv3"] <- rnorm(200, 70, 20)
#othvars[sample(1:200, size=5), "iv1"] <- NA
#othvars[sample(1:200, size=5), "iv2"] <- NA
#othvars[sample(1:200, size=5), "cv1"] <- NA
#othvars[sample(1:200, size=5), "cv2"] <- NA
#othvars[sample(1:200, size=5), "cv3"] <- NA

#testdata <- join(testdata, othvars, by = "id")
#testdata <- dcast(testdata, id + iv1 + iv2 + cv1 + cv2 + cv3 ~ variable + time)

# Test cases #


#rm.ssp(df=testdata,
#       .factor1=list(A,locationB,locationC,locationD,locationE),
#       depvars=list(dv1_1,dv1_2,dv1_3,dv1_4,dv1_5),
#       indepvars=list(),
#       indepvars.int=list(),
#       covars=list(),
#       covars.int=list()) -> a

#rm.ssp(df=testdata,
##       .factor1=list(locationA,locationB,locationC,locationD,locationE),
#       depvars=list(dv1_1,dv1_2,dv1_3,dv1_4,dv1_5),
#       indepvars=list(iv1),
#       indepvars.int=list(iv1:iv2),
#       covars=list(),
#       covars.int=list()) -> a

#rm.ssp(df=testdata,
#       .factor1=list(locationA,locationB,locationC,locationD,locationE),
#       depvars=list(dv1_1,dv1_2,dv1_3,dv1_4,dv1_5),
#       indepvars=list(iv1),
#       indepvars.int=list(iv1:iv2),
#       covars=list(cv1),
#       covars.int=list(cv2:cv3)) -> a

#rm.ssp(df=testdata,
#       .factor1=list(locationA,locationB,locationC,locationD,locationE),
#       depvars=list(dv1_1,dv1_2,dv1_3,dv1_4,dv1_5),
#       indepvars=list(iv1,iv2),
#       indepvars.int=list(iv1:iv2),
#       covars=list(),
#       covars.int=list()) -> a

#rm.ssp(df=testdata,
#       .factor1=list(locationA,locationB,locationC,locationD,locationE),
#       depvars=list(dv1_1,dv1_2,dv1_3,dv1_4,dv1_5),
#       indepvars=list(),
#       indepvars.int=list(),
#       covars=list(cv1),
#      covars.int=list(cv1:cv2)) -> a

#rm.ssp(df=testdata,
#       .factor1=list(locationA,locationB,locationC,locationD,locationE),
#       depvars=list(dv1_1,dv1_2,dv1_3,dv1_4,dv1_5),
#       indepvars=list(),
#       indepvars.int=list(),
#       covars=list(cv1,cv2),
#       covars.int=list(cv1:cv2)) -> a

#rm.ssp(df=testdata,
#            .factor1=list(locationA,locationB,locationC,locationD,locationE),
#            depvars=list(dv1_1,dv1_2,dv1_3,dv1_4,dv1_5),
#            indepvars=list(iv1,iv2),
#            indepvars.int=list(iv1:iv2),
#            covars=list(cv1,cv2,cv3),
#            covars.int=list(cv1:cv3)) -> a

#rm.ssp(df=testdata,
#       .factor1=list(locationA,locationB,locationC,locationD,locationE),
#       .factor2=list(pretest,posttest,followup),
#       depvars=list(dv1_1,dv1_2,dv1_3,dv1_4,dv1_5,
#                    dv2_1,dv2_2,dv2_3,dv2_4,dv2_5,
#                    dv3_1,dv3_2,dv3_3,dv3_4,dv3_5),
#       indepvars=list(iv1,iv2),
#       indepvars.int=list(iv1:iv2),
#       covars=list(cv1,cv2,cv3),
#       covars.int=list(cv1:cv3)) -> a

#rm.ssp(df=testdata,
#       .factor1=list(locationA,locationB,locationC,locationD,locationE),
#       .factor2=list(pretest,posttest,followup),
#       depvars=list(dv1_1,dv1_2,dv1_3,dv1_4,dv1_5,
#                    dv2_1,dv2_2,dv2_3,dv2_4,dv2_5,
#                    dv3_1,dv3_2,dv3_3,dv3_4,dv3_5),
#       indepvars=list(),
#       indepvars.int=list(),
#       covars=list(cv1,cv2,cv3),
#       covars.int=list(cv1:cv3)) -> a

#rm.ssp(df=testdata,
#       .factor1=list(locationA,locationB,locationC,locationD,locationE),
#       .factor2=list(pretest,posttest,followup),
#       depvars=list(dv1_1,dv1_2,dv1_3,dv1_4,dv1_5,
#                    dv2_1,dv2_2,dv2_3,dv2_4,dv2_5,
#                    dv3_1,dv3_2,dv3_3,dv3_4,dv3_5),
#       indepvars=list(),
#       indepvars.int=list(),
#       covars=list(),
#       covars.int=list()) -> a

#rm.ssp(df=testdata,
#       .factor1=list(locationA,locationB,locationC,locationD,locationE),
#       .factor2=list(pretest,posttest,followup),
#       .factor.int=list(factor1:iv1,factor2:iv2),
#       depvars=list(dv1_1,dv1_2,dv1_3,dv1_4,dv1_5,
#                    dv2_1,dv2_2,dv2_3,dv2_4,dv2_5,
#                    dv3_1,dv3_2,dv3_3,dv3_4,dv3_5),
#       indepvars=list(iv1,iv2),
#       indepvars.int=list(iv1:iv2),
#       covars=list(cv1,cv2),
#       covars.int=list(cv1:cv2,iv1:cv2)) -> a





#data <- read.csv(file=file.choose())
#RM ANCOVA. 3 DVs, 1 IV, 1 CV
#b <- rm.ssp(df=data,
#            .factor1=list(Normal.Dependent, Not.Normal.Dependent, Normal),
#            depvars=list(Normal.Dependent, Not.Normal.Dependent, Normal),
#            indepvars=list(Group.3),
#            covars=list(IVModerate))
#b

#RM ANCOVA, 3 DVs, 0 IV, 1 CV
#b <- rm.ssp(df=data,
#            .factor1=list(Normal.Dependent, Not.Normal.Dependent, Normal),
#            depvars=list(Normal.Dependent, Not.Normal.Dependent, Normal),
#            indepvars=list(),
#            covars=list(IVModerate))
#b


#testdata <- melt(testdata, id = c("id", "time"))
#testdata <- testdata[-sample(nrow(testdata), size = 50, replace = FALSE),]