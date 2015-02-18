#' MANOVA and MANCOVA
#' Required Libraries : car, plyr, MVN, biotools, lsmeans, data.table 
#' 
#' @Perform  MANOVA or MANCOVA 
#' @param df - data frame that contains all the required variables
#' @param depvars - dependent variables, can be repeated measures (within subject variables)
#' @param indepvars - independent variables or factors (between subject variables)
#' @param indepvars.int - interaction terms of indepvars 
#' @param covars - covariates that reduce the error variance, continuous
#' @param covars.int - interaction terms of covars or covars:indepvars  
#' @return output - parameters and statistics

#options(digits=4)

#library(car) # Anova #
#library(plyr)
#library(MVN) # Mardia #
#library(biotools) # Box's M #
#library(lattice) # splom #
#library(lsmeans) # marginal means #
#library(data.table)
#library(jsonlite)


##################################################################################################################################
################################################## Start of Main Function ########################################################
##################################################################################################################################

mancova.ssp <- function(df,
                        depvars = list(),
                        indepvars = list(),
                        indepvars.int = list(),
                        covars = list(),
                        covars.int = list()) {
        
        ###############################
        # Check for duplicate columns #
        ###############################
        
        # Convert variable lists into vectors of variable names for assumption test #
        # Make sure every variable specified in the parameters are in df.complete #
        depvars <- paste(substitute(depvars))[-1]
        indepvars <- paste(substitute(indepvars))[-1]
        indepvars.int <- paste(substitute(indepvars.int))[-1]
        covars <- paste(substitute(covars))[-1]
        covars.int <- paste(substitute(covars.int))[-1]
        
        # 1) Use complete obs based on variable specification in the parameters #
        # 2) Indicate which columns are duplicate and to be removed #
        # 3) Continue running the codes with unique columns #
        
        if (length(indepvars.int)>0)
        {all.indepvars <- unlist(strsplit(indepvars.int,split=":"))
         all.indepvars <- unique(c(indepvars, all.indepvars))
        } else {all.indepvars <- indepvars} 
        if (length(covars.int)>0)
        {all.covars <- unlist(strsplit(covars.int,split=":"))
         all.covars <- unique(c(covars, all.covars))
        } else {all.covars <- covars}   
        
        df.complete <- df[, c(depvars, all.indepvars, all.covars)]
        df.complete <- df.complete[complete.cases(df.complete),]
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
                    xlab = "",
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
        #boxm.test.results <- list()
        #for (iv in indepvars) {
        #        boxm.test.results[[iv]] <- boxM(df.complete[,depvars], df.complete[,iv])
        #}
        #boxm.final.results <- data.frame(matrix(unlist(lapply(boxm.test.results, function(x) c(x$statistic, x$parameter, x$p.value))),ncol=3,byrow=TRUE),stringsAsFactors=FALSE)  
        #names(boxm.final.results) <- c("Chi-square Statistics","df","p.value")
        #boxm.final.results[,"depvar"] <- paste(depvars, collapse = ",")
        #boxm.final.results[,"indepvar"] <- indepvars
        
        
        # 5) Homogeneity of regression slopes (between DVs and covariates), if covariates are provided #
        if (length(covars)>0) {
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
        
        
        # 7) IV frequency #
        iv.freq <- lapply(as.data.frame(df.complete[,indepvars]), function(x) as.data.frame(table(x)))
        
        ####################
        # MANOVA / MANCOVA #
        ####################
        
        # Process divided into two parts #
        # 1) Run MANOVA, determine marginal means, multivariate and univariate outcomes #
        # 2) If covariates are provided, run MANCOVA to get the above outcomes #
        
        # Variables to be inserted into the lm formula # 
        depvarlist <- paste(depvars, collapse = ",")
        indepvarlist <- paste(c(indepvars, indepvars.int), collapse = "+")
        
        ########################################
        # Analyses without covariates (MANOVA) #
        ########################################
        
        manova.results <- eval(parse(text=paste0("lm(cbind(", depvarlist, ") ~ ", indepvarlist, ", data = df.complete)"))) 
        
        # Output Multivariate test results #
        outtests <- car:::print.Anova.mlm
        # allow the function to return the results and disable print
        body(outtests)[[16]] <- quote(invisible(tests))
        body(outtests)[[15]] <- NULL
        # Run the Anova over all tests  
        mtests <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
        manova.mtests.results <- lapply(mtests, function(i) outtests(Anova(manova.results, test.statistic=i,type=3))[-1,])
        # Stack all results into one data frame #
        terms <- row.names(manova.mtests.results[[1]])
        multivariate.test <- rep(mtests, lapply(manova.mtests.results, function(x) dim(x)[1]))
        manova.mtests.results <- as.data.frame(do.call("rbind", manova.mtests.results))  
        manova.mtests.results[,"term"] <- terms
        manova.mtests.results[,"Multivariate test"] <- multivariate.test  
        row.names(manova.mtests.results) <- NULL
        
        # Univariate models to calculate marginal means #
        manova.models <- lapply(depvars, function(x) as.formula(paste(x, " ~ ", indepvarlist)))
        manova.models <- lapply(manova.models, function(x) lm(x, data = df.complete))
        
        # Univariate analysis #
        manova.univariate <- lapply(manova.models, function(x) Anova(x, type = 3))
        response <- rep(depvars, lapply(manova.univariate, function(x) dim(x)[1]))
        manova.univariate <- as.data.frame(do.call("rbind", manova.univariate))
        manova.univariate[,"Mean Sq"] <- manova.univariate[,"Sum Sq"]/manova.univariate[,"Df"]
        manova.univariate[,"depvar"] <- response
        manova.univariate[,"indepvar"] <- c("Intercept", c(indepvars, indepvars.int), "Residuals")
        row.names(manova.univariate) <- NULL
        
        # Marginal means #
        # Default contrast - pairwise #
        # Grouping variables (indepvars) have to be factors #
        # Sample size of each group must be at least 3 #
        manova.lsmeans <- unlist(lapply(manova.models, function(x) sapply(indepvars, function(y) lsmeans(x,y))))
        
        # Post hoc analysis - Pairwise Comparison using Tukey HSD on Marginal means #
        manova.posthoc <- lapply(manova.lsmeans, function(x) pairs(x)) 
        # Retrieve sig level #
        manova.posthoc.sig <- lapply(manova.posthoc, function(x) summary(x)) 
        # Retrieve confidence interval
        manova.posthoc.ci <- lapply(manova.lsmeans, function(x) confint(pairs(x)))
        
        # Combine outputs #
        grpname <- rep(names(manova.posthoc.sig),lapply(manova.posthoc.sig, function(x) dim(x)[1])) 
        manova.posthoc.sig <- as.data.frame(rbindlist(manova.posthoc.sig))
        manova.posthoc.ci <- as.data.frame(rbindlist(manova.posthoc.ci)) 
        manova.posthoc.final <- cbind(manova.posthoc.sig, manova.posthoc.ci[, c("lower.CL","upper.CL")])
        manova.posthoc.final[,"depvar"] <- rep(depvars, each = nrow(manova.posthoc.final)/length(depvars))
        manova.posthoc.final[,"indepvar"] <- grpname
        
        
        # Output Marginal means to a data frame #
        # Change variable names and stack all data frames together #
        manova.lsmeans <- lapply(manova.lsmeans, function(x) summary(x))
        grpname <- rep(names(manova.lsmeans),lapply(manova.lsmeans, function(x) dim(x)[1]))
        manova.lsmeans <- as.data.frame(rbindlist(manova.lsmeans))
        manova.lsmeans[,"depvar"] <- rep(depvars, each = nrow(manova.lsmeans)/length(depvars))
        manova.lsmeans[,"indepvar"] <- grpname
        names(manova.lsmeans)[1] <- "indepvarvalue"
        
        
        ######################################
        # Analyses with covariates (MANCOVA) #
        ######################################
        
        if (length(covars)==0 & length(covars.int)==0) {
                mancova.mtests.results <- mancova.univariate <- mancova.lsmeans <- mancova.posthoc.final <- NULL
        } else { 
                covarlist <- paste(c(covars, covars.int), collapse="+")
                mancova.results <- eval(parse(text = paste0("lm(cbind(", depvarlist, ") ~ ", indepvarlist, " + ", covarlist, ", data = df.complete)")))
                
                # Output Multivariate test results #
                mancova.mtests.results <- lapply(mtests, function(i) outtests(Anova(mancova.results, test.statistic=i,type=3))[-1,])
                # Stack all results into one data frame #
                terms <- row.names(mancova.mtests.results[[1]])
                multivariate.test <- rep(mtests, lapply(mancova.mtests.results, function(x) dim(x)[1]))
                mancova.mtests.results <- as.data.frame(do.call("rbind", mancova.mtests.results))  
                mancova.mtests.results[,"term"] <- terms
                mancova.mtests.results[,"Multivariate test"] <- multivariate.test  
                row.names(mancova.mtests.results) <- NULL
                
                # Univariate models to calculate marginal means #
                mancova.models <- lapply(depvars, function(x) as.formula(paste(x," ~ ", indepvarlist, "+", covarlist)))
                mancova.models <- lapply(mancova.models, function(x) lm(x, data = df.complete))
                
                # Univariate analysis #
                mancova.univariate <- lapply(mancova.models, function(x) Anova(x, type = 3))
                response <- rep(depvars, lapply(mancova.univariate, function(x) dim(x)[1]))
                mancova.univariate <- as.data.frame(do.call("rbind", mancova.univariate))
                mancova.univariate[,"Mean Sq"] <- mancova.univariate[,"Sum Sq"]/mancova.univariate[,"Df"]
                mancova.univariate[,"depvar"] <- response
                mancova.univariate[,"indepvar"] <- c("Intercept", c(indepvars, indepvars.int, covars, covars.int), "Residuals")
                row.names(mancova.univariate) <- NULL
                
                # Marginal means #
                # Default contrast - pairwise #
                # Grouping variables (indepvars) have to be factors #
                # Sample size of each group must be at least 3 #
                mancova.lsmeans <- unlist(lapply(mancova.models, function(x) sapply(indepvars, function(y) lsmeans(x,y))))
                
                # Post hoc analysis - Pairwise Comparison using Tukey HSD on Marginal means #
                mancova.posthoc <- lapply(mancova.lsmeans, function(x) pairs(x))
                # Retrieve sig level #
                mancova.posthoc.sig <- lapply(mancova.posthoc, function(x) summary(x))
                # Retrieve confidence interval
                mancova.posthoc.ci <- lapply(mancova.lsmeans, function(x) confint(pairs(x)))
                
                # Combine outputs #
                grpname <- rep(names(mancova.posthoc.sig),lapply(mancova.posthoc.sig, function(x) dim(x)[1]))
                mancova.posthoc.sig <- as.data.frame(rbindlist(mancova.posthoc.sig))
                mancova.posthoc.ci <- as.data.frame(rbindlist(mancova.posthoc.ci))
                mancova.posthoc.final <- cbind(mancova.posthoc.sig, mancova.posthoc.ci[, c("lower.CL","upper.CL")])
                mancova.posthoc.final[,"depvar"] <- rep(depvars, each = nrow(mancova.posthoc.final)/length(depvars))
                mancova.posthoc.final[,"indepvar"] <- grpname
                
                # Output Marginal means to a data frame #
                # Change variable names and stack all data frames together #
                mancova.lsmeans <- lapply(mancova.lsmeans, function(x) summary(x))
                grpname <- rep(names(mancova.lsmeans),lapply(mancova.lsmeans, function(x) dim(x)[1]))
                mancova.lsmeans <- as.data.frame(rbindlist(mancova.lsmeans))
                mancova.lsmeans[,"depvar"] <- rep(depvars, each = nrow(mancova.lsmeans)/length(depvars))
                mancova.lsmeans[,"indepvar"] <- grpname
                names(mancova.lsmeans)[1] <- "indepvarvalue"   
                
        }
        
        
        # Final Output #  
        return(list(uvn = uvn.final.results, mvn = mvn.final.results, 
                    bartlett = bartlett.final.results, 
                    levene = levene.final.results, 
                    slopes = equal.slopes, corr = dvcorr, n = iv.freq, 
                    manova = manova.mtests.results, anova = manova.univariate, adj.manova = manova.lsmeans, manova.post = manova.posthoc.final,
                    mancova = mancova.mtests.results, ancova = mancova.univariate, adj.mancova = mancova.lsmeans, mancova.post = mancova.posthoc.final)) 
}
  

##################################################################################################################################
################################################## End of Main Function ##########################################################
##################################################################################################################################

# Testing #

#testdata <- read.csv("testdata.csv")
# Converted IV (grouping variables) to character #
#testdata[,"TypicalLikert1"] <- factor(testdata[,"TypicalLikert1"])
#testdata[,"TypicalLikert2"] <- factor(testdata[,"TypicalLikert2"])
#testdata[,"TypicalLikert3"] <- factor(testdata[,"TypicalLikert3"])
#testdata[,"TypicalLikert9"] <- factor(sample(c(1:7),size=nrow(testdata),replace=TRUE))

# Testing without missing data #
#mancova.ssp(df=testdata,
 #           depvars=list(NormalDependent,Normal,Paired4),
  #          indepvars=list(TypicalLikert1,BigGroup,TypicalLikert9),
   #         indepvars.int=list(ChiGroup1:BigGroup),
    #        covars=list(IVModerate,Paired5,Paired6),
     #       covars.int=list(TypicalLikert1:IVModerate,TypicalLikert9:Paired6,Paired5:Paired6)) -> a


#mancova.ssp(df=testdata,
 #           depvars=list(NormalDependent,Normal,Paired4),
  #          indepvars=list(TypicalLikert1,BigGroup,TypicalLikert9),
   #         covars=list(IVModerate,Paired5,Paired6),
    #        covars.int=list(TypicalLikert1:IVModerate,TypicalLikert9:Paired6,Paired5:Paired6)) -> a

#mancova.ssp(df=testdata,
 #           depvars=list(NormalDependent,Normal,Paired4),
  #          indepvars=list(TypicalLikert1,BigGroup,TypicalLikert9),
   #         covars=list(IVModerate,Paired5,Paired6)) -> a

     
# Reminders #
# 0) Users should be reminded to remove outliers as this procedure is very sensitive to outliers
# 1) Independent vars have to be factors
# 2) Inverse of Model Matrix doesn't exist (BigGroup, TypicalLikert3) interacts with other independent variables - Some zeros in the cross tab - Design Matrix columns = 0#
#    -  Error in solve.default(crossprod(model.matrix(mod))) : 
#       Lapack routine dgesv: system is exactly singular: U[36,36] = 0 -> Provide a warning of linear depedence between variables
# 3) How to input interaction terms - need a discussion on interface 
# 4) Marginal Means are used instead of descriptive means with pairwise contrast as default for post hoc comparison using TukeyHSD
# 5) descriptive stats of DV by each IV level - existing functions? #
# 6) Only complete cases will be used
# 7) All F-tests use Type III Sum of squares


#data <- read.csv(file=file.choose())

#a <- mancova.ssp(df=data, depvars=list(Paired.1, Paired.2.Different, Normal), indepvars=list(Big.Group), covars=list(IVModerate))
#a
#fromJSON(toJSON(a))

#b <- mancova.ssp(df=data, depvars=list(Paired.1, Paired.2.Different, Normal), indepvars=list(Different.Group.2, Group.3))
#b
#fromJSON(toJSON(b))

