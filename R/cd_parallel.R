# Reproducing the observed correlation matrix rather than generating random data #

########################################################################################################################
GenData <- function(Supplied.Data, N.Factors, N, Max.Trials = 5, Initial.Multiplier = 1, Cor.Type)
{
  # Steps refer to description in the following article:
  # Ruscio, J., & Kaczetow, W. (2008). Simulating multivariate nonnormal data using an iterative algorithm. 
  # Multivariate Behavioral Research, 43(3), 355-381.
  
  # Initialize variables and (if applicable) set random number seed (step 1) -------------------------------------
  
  set.seed(1)
  k <- dim(Supplied.Data)[2]
  Data <- matrix(0, nrow = N, ncol = k)            # Matrix to store the simulated data
  Distributions <- matrix(0, nrow = N, ncol = k)   # Matrix to store each variable's score distribution
  Iteration <- 0                                   # Iteration counter
  Best.RMSR <- 1                                   # Lowest RMSR correlation
  Trials.Without.Improvement <- 0                  # Trial counter
  
  # Generate distribution for each variable (step 2) -------------------------------------------------------------
  
  for (i in 1:k)
    Distributions[,i] <- sort(sample(Supplied.Data[,i], size = N, replace = T))
  
  # Calculate and store a copy of the target correlation matrix (step 3) -----------------------------------------
  
  Target.Corr <- cor(Supplied.Data, method = Cor.Type)
  Intermediate.Corr <- Target.Corr
  
  # Generate random normal data for shared and unique components, initialize factor loadings (steps 5, 6) --------
  
  Shared.Comp <- matrix(rnorm(N * N.Factors, 0, 1), nrow = N, ncol = N.Factors)
  Unique.Comp <- matrix(rnorm(N * k, 0, 1), nrow = N, ncol = k)
  Shared.Load <- matrix(0, nrow = k, ncol = N.Factors)
  Unique.Load <- matrix(0, nrow = k, ncol = 1)
  
  # Begin loop that ends when specified number of iterations pass without improvement in RMSR correlation --------
  
  while (Trials.Without.Improvement < Max.Trials)
  {
    Iteration <- Iteration + 1
    
    # Calculate factor loadings and apply to reproduce desired correlations (steps 7, 8) ---------------------------
    
    Fact.Anal <- Factor.Analysis(Intermediate.Corr, Corr.Matrix = TRUE, N.Factors = N.Factors, Cor.Type = Cor.Type)
    if (N.Factors == 1) Shared.Load[,1] <- Fact.Anal$loadings
    else 
      for (i in 1:N.Factors)
        Shared.Load[,i] <- Fact.Anal$loadings[,i]
    Shared.Load[Shared.Load > 1] <- 1
    Shared.Load[Shared.Load < -1] <- -1
    if (Shared.Load[1,1] < 0) Shared.Load <- Shared.Load * -1
    for (i in 1:k)
      if (sum(Shared.Load[i,] * Shared.Load[i,]) < 1) Unique.Load[i,1] <- 
      (1 - sum(Shared.Load[i,] * Shared.Load[i,]))
    else Unique.Load[i,1] <- 0
    Unique.Load <- sqrt(Unique.Load)
    for (i in 1:k)
      Data[,i] <- (Shared.Comp %*% t(Shared.Load))[,i] + Unique.Comp[,i] * Unique.Load[i,1]
    
    # Replace normal with nonnormal distributions (step 9) ---------------------------------------------------------
    
    for (i in 1:k)
    {
      Data <- Data[sort.list(Data[,i]),]
      Data[,i] <- Distributions[,i]
    }
    
    # Calculate RMSR correlation, compare to lowest value, take appropriate action (steps 10, 11, 12) --------------
    
    Reproduced.Corr <- cor(Data, method = Cor.Type)
    Residual.Corr <- Target.Corr - Reproduced.Corr
    RMSR <- sqrt(sum(Residual.Corr[lower.tri(Residual.Corr)] * Residual.Corr[lower.tri(Residual.Corr)]) / 
                   (.5 * (k * k - k)))
    if (RMSR < Best.RMSR)
    {
      Best.RMSR <- RMSR
      Best.Corr <- Intermediate.Corr
      Best.Res <- Residual.Corr
      Intermediate.Corr <- Intermediate.Corr + Initial.Multiplier * Residual.Corr
      Trials.Without.Improvement <- 0
    }
    else 
    {
      Trials.Without.Improvement <- Trials.Without.Improvement + 1
      Current.Multiplier <- Initial.Multiplier * .5 ^ Trials.Without.Improvement
      Intermediate.Corr <- Best.Corr + Current.Multiplier * Best.Res
    }
  }
  
  # Construct the data set with the lowest RMSR correlation (step 13) --------------------------------------------
  
  Fact.Anal <- Factor.Analysis(Best.Corr, Corr.Matrix = TRUE, N.Factors = N.Factors, Cor.Type = Cor.Type)
  if (N.Factors == 1) Shared.Load[,1] <- Fact.Anal$loadings
  else
    for (i in 1:N.Factors)
      Shared.Load[,i] <- Fact.Anal$loadings[,i]
  Shared.Load[Shared.Load > 1] <- 1
  Shared.Load[Shared.Load < -1] <- -1
  if (Shared.Load[1,1] < 0) Shared.Load <- Shared.Load * -1
  for (i in 1:k)
    if (sum(Shared.Load[i,] * Shared.Load[i,]) < 1) Unique.Load[i,1] <-
    (1 - sum(Shared.Load[i,] * Shared.Load[i,]))
  else Unique.Load[i,1] <- 0
  Unique.Load <- sqrt(Unique.Load)
  for (i in 1:k)
    Data[,i] <- (Shared.Comp %*% t(Shared.Load))[,i] + Unique.Comp[,i] * Unique.Load[i,1]
  Data <- apply(Data, 2, scale) # standardizes each variable in the matrix
  for (i in 1:k)
  {
    Data <- Data[sort.list(Data[,i]),]
    Data[,i] <- Distributions[,i]
  }
  
  # Return the simulated data set (step 14) ----------------------------------------------------------------------
  
  return(Data)
}

####################################################################################################################
Factor.Analysis <- function(Data, Corr.Matrix = FALSE, Max.Iter = 50, N.Factors = 0, Cor.Type)
{
  Data <- as.matrix(Data)
  k <- dim(Data)[2]
  if (N.Factors == 0)
  {
    N.Factors <- k
    Determine <- T
  }
  else Determine <- F
  if (!Corr.Matrix) Cor.Matrix <- cor(Data, method = Cor.Type)
  else Cor.Matrix <- Data
  Criterion <- .001
  Old.H2 <- rep(99, k)
  H2 <- rep(0, k)
  Change <- 1
  Iter <- 0
  Factor.Loadings <- matrix(nrow = k, ncol = N.Factors)
  while ((Change >= Criterion) & (Iter < Max.Iter))
  {
    Iter <- Iter + 1
    Eig <- eigen(Cor.Matrix)
    L <- sqrt(Eig$values[1:N.Factors])
    for (i in 1:N.Factors)
      Factor.Loadings[,i] <- Eig$vectors[,i] * L[i]
    for (i in 1:k)
      H2[i] <- sum(Factor.Loadings[i,] * Factor.Loadings[i,])
    Change <- max(abs(Old.H2 - H2))
    Old.H2 <- H2
    diag(Cor.Matrix) <- H2
  }
  if (Determine) N.Factors <- sum(Eig$values > 1)
  return(list(loadings = Factor.Loadings[,1:N.Factors], factors = N.Factors))
}

##########################################################################################################################
################################################### Comparison Data ######################################################
##########################################################################################################################
fa.cd <- function(Data, F.Max, N.Pop = 10000, N.Samples = 500, Alpha = .30, Graph = F, Spearman = F)
{
  # Data = N (sample size) x k (number of variables) data matrix
  # F.Max = largest number of factors to consider
  # N.Pop = size of finite populations of comparison data (default = 10,000 cases)
  # N.Samples = number of samples drawn from each population (default = 500)
  # Alpha = alpha level when testing statistical significance of improvement with add'l factor (default = .30) 
  # Graph = whether to plot the fit of eigenvalues to those for comparison data (default = F)
  # Spearman = whether to use Spearman rank-order correlations rather than Pearson correlations (default = F)
  # Use complete.obs because of resampling #
  
  N <- dim(Data)[1]
  k <- dim(Data)[2]
  if (Spearman) Cor.Type <- "spearman" else Cor.Type <- "pearson"
  cor.Data <- cor(x=Data, method = Cor.Type)
  Eigs.Data <- eigen(cor.Data)$values
  RMSR.Eigs <- matrix(0, nrow = N.Samples, ncol = F.Max)
  Sig <- T
  F.CD <- 1
  while ((F.CD <= F.Max) & (Sig))
  {
    Pop <- GenData(Data, N.Factors = F.CD, N = N.Pop, Cor.Type = Cor.Type)
    for (j in 1:N.Samples)
    {
      Samp <- Pop[sample(1:N.Pop, size = N, replace = T),]
      cor.Samp <- cor(Samp, method = Cor.Type)
      Eigs.Samp <- eigen(cor.Samp)$values
      RMSR.Eigs[j,F.CD] <- sqrt(sum((Eigs.Samp - Eigs.Data) * (Eigs.Samp - Eigs.Data)) / k)
    }
    if (F.CD > 1) Sig <- (wilcox.test(RMSR.Eigs[,F.CD], RMSR.Eigs[,(F.CD - 1)], "less")$p.value < Alpha)
    if (Sig) F.CD <- F.CD + 1
  }
  if (Graph)
  {
    if (Sig) x.max <- F.CD - 1
    else x.max <- F.CD
    ys <- apply(RMSR.Eigs[,1:x.max], 2, mean)
    plot(x = 1:x.max, y = ys, ylim = c(0, max(ys)), xlab = "Factor Number", ylab = "RMSR Eigenvalue", type = "b", 
         main = "Fit to Comparison Data")
    abline(v = F.CD - 1, lty = 3)
  }
  nfact <- F.CD - 1
  return(nfact)
}

####################################################################################################################
################################ Modified parallel analysis ########################################################
####################################################################################################################

fa.parallel1 <- function (x, n.obs = NULL, fm = "minres", fa = "both", main = "Parallel Analysis Scree Plots", 
                          n.iter = 20, error.bars = FALSE, SMC = FALSE, ylabel = NULL, 
                          show.legend = TRUE, sim = TRUE, use) 
{
  cl <- match.call()
  if (!require(parallel)) {
    message("The parallel package needs to be installed to run mclapply")
  }
  ci <- 1.96
  arrow.len <- 0.05
  nsub <- dim(x)[1]
  nvariables <- dim(x)[2]
  if ((nsub == nvariables) && !sim) {
    warning("You specified a correlation matrix, but asked to just resample (sim was set to FALSE).  This is impossible, so sim is set to TRUE")
    sim <- TRUE
  }
  if (!is.null(n.obs)) {
    nsub <- n.obs
    rx <- x
    if (dim(x)[1] != dim(x)[2]) {
      warning("You specified the number of subjects, implying a correlation matrix, but do not have a correlation matrix, correlations found ")
      rx <- cor(x, use = use)
      if (!sim) {
        warning("You specified a correlation matrix, but asked to just resample (sim was set to FALSE).  This is impossible, so sim is set to TRUE")
        sim <- TRUE
      }
    }
  }
  else {
    if (nsub == nvariables) {
      warning("It seems as if you are using a correlation matrix, but have not specified the number of cases. The number of subjects is arbitrarily set to be 100  ")
      rx <- x
      nsub = 100
      n.obs = 100
    }
    else {
      rx <- cor(x, use = use)
    }
  }
  valuesx <- eigen(rx)$values
  if (SMC) {
    diag(rx) <- smc(rx)
    fa.valuesx <- eigen(rx)$values
  }
  else {
    fa.valuesx <- fa(rx, fm = fm, warnings = FALSE)$values
  }
  temp <- list(samp = vector("list", n.iter), samp.fa = vector("list", 
                                                               n.iter), sim = vector("list", n.iter), sim.fa = vector("list", 
                                                                                                                      n.iter))
  templist <- mclapply(1:n.iter, function(XX) {
    if (is.null(n.obs)) {
      bad <- TRUE
      while (bad) {
        sampledata <- matrix(apply(x, 2, function(y) sample(y, 
                                                            nsub, replace = TRUE)), ncol = nvariables)
        C <- cor(sampledata, use = use)
        bad <- any(is.na(C))
      }
      values.samp <- eigen(C)$values
      temp[["samp"]] <- values.samp
      if (fa != "pc") {
        if (SMC) {
          sampler <- C
          diag(sampler) <- smc(sampler)
          temp[["samp.fa"]] <- eigen(sampler)$values
        }
        else {
          temp[["samp.fa"]] <- fa(C, fm = fm, SMC = FALSE, 
                                  warnings = FALSE)$values
        }
      }
    }
    if (sim) {
      simdata = matrix(rnorm(nsub * nvariables), nrow = nsub, 
                       ncol = nvariables)
      sim.cor <- cor(simdata)
      temp[["sim"]] <- eigen(sim.cor)$values
      if (fa != "pc") {
        if (SMC) {
          diag(sim.cor) <- smc(sim.cor)
          temp[["sim.fa"]] <- eigen(sim.cor)$values
        }
        else {
          fa.values.sim <- fa(sim.cor, fm = fm, SMC = FALSE, 
                              warnings = FALSE)$values
          temp[["sim.fa"]] <- fa.values.sim
        }
      }
    }
    replicates <- list(samp = temp[["samp"]], samp.fa = temp[["samp.fa"]], 
                       sim = temp[["sim"]], sim.fa = temp[["sim.fa"]])
  })
  if (is.null(ylabel)) {
    if (fa != "pc") {
      ylabel <- "Eigenvalues of Principal Components and Factor Analysis"
    }
    else {
      ylabel <- "Eigen values of Principal Components"
    }
  }
  values <- t(matrix(unlist(templist), ncol = n.iter))
  values.sim.mean = colMeans(values, na.rm = TRUE)
  values.sim.se <- apply(values, 2, sd, na.rm = TRUE)/sqrt(n.iter)
  ymax <- max(valuesx, values.sim.mean)
  if (sim) {
    sim.pc <- values.sim.mean[1:nvariables]
    sim.fa <- values.sim.mean[(nvariables + 1):(2 * nvariables)]
    sim.pcr <- values.sim.mean[(2 * nvariables + 1):(3 * 
                                                       nvariables)]
    sim.far <- values.sim.mean[(3 * nvariables + 1):(4 * 
                                                       nvariables)]
  }
  else {
    sim.pcr <- values.sim.mean[1:nvariables]
    sim.far <- values.sim.mean[(nvariables + 1):(2 * nvariables)]
    sim.fa <- NA
    sim.pc <- NA
  }
  if (fa != "fa") {
    if (fa == "both") {
      plot(valuesx, type = "b", main = main, ylab = ylabel, 
           ylim = c(0, ymax), xlab = "Factor/Component Number", 
           pch = 4, col = "black")
    }
    else {
      plot(valuesx, type = "b", main = main, ylab = ylabel, 
           ylim = c(0, ymax), xlab = "Component Number", 
           pch = 4, col = "black")
    }
    if (sim) 
      points(sim.pc, type = "l", lty = "dotted", pch = 4, 
             col = "black")
    if (error.bars) {
      sim.se.pc <- values.sim.se[1:nvariables]
      sim.se.fa <- values.sim.se[(nvariables + 1):(2 * 
                                                     nvariables)]
      sim.pcr.se <- values.sim.se[(2 * nvariables + 1):(3 * 
                                                          nvariables)]
      sim.se.fa <- values.sim.se[(3 * nvariables + 1):(4 * 
                                                         nvariables)]
      for (i in 1:length(sim.pc)) {
        ycen <- sim.pc[i]
        yse <- sim.se.pc[i]
        arrows(i, ycen - ci * yse, i, ycen + ci * yse, 
               length = arrow.len, angle = 90, code = 3, col = par("fg"), 
               lty = NULL, lwd = par("lwd"), xpd = NULL)
      }
    }
    if (is.null(n.obs)) {
      points(sim.pcr, type = "l", lty = "dashed", pch = 4, 
             col = "black")
      if (error.bars) {
        for (i in 1:length(sim.pcr)) {
          ycen <- sim.pcr[i]
          yse <- sim.pcr.se[i]
          arrows(i, ycen - ci * yse, i, ycen + ci * yse, 
                 length = arrow.len, angle = 90, code = 3, 
                 col = par("fg"), lty = NULL, lwd = par("lwd"), 
                 xpd = NULL)
        }
      }
    }
  }
  if (fa != "pc") {
    if (fa == "fa") {
      ylabel <- "Eigen Values of Principal Factors"
      plot(fa.valuesx, type = "b", main = main, ylab = ylabel, 
           ylim = c(0, ymax), xlab = "Factor Number", pch = 4, 
           col = "black")
    }
    if (sim) 
      points(sim.fa, type = "l", lty = "dotted", pch = 2, 
             col = "black")
    if (error.bars) {
      for (i in 1:length(sim.fa)) {
        ycen <- sim.fa[i]
        yse <- sim.se.fa[i]
        arrows(i, ycen - ci * yse, i, ycen + ci * yse, 
               length = arrow.len, angle = 90, code = 3, col = par("fg"), 
               lty = NULL, lwd = par("lwd"), xpd = NULL)
      }
    }
    if (is.null(n.obs)) {
      points(sim.far, type = "l", lty = "dashed", pch = 2, 
             col = "black")
      if (error.bars) {
        for (i in 1:length(sim.far)) {
          ycen <- sim.far[i]
          yse <- sim.se.fa[i]
          arrows(i, ycen - ci * yse, i, ycen + ci * yse, 
                 length = arrow.len, angle = 90, code = 3, 
                 col = par("fg"), lty = NULL, lwd = par("lwd"), 
                 xpd = NULL)
        }
      }
    }
    if (fa != "fa") 
      points(fa.valuesx, type = "b", lty = "solid", pch = 2, 
             col = "black")
    if (sim) 
      points(sim.fa, type = "l", lty = "dotted", pch = 2, 
             col = "black")
    if (is.null(n.obs)) {
      points(sim.far, type = "l", lty = "dashed", pch = 2, 
             col = "black")
    }
  }
  if (show.legend) {
    if (is.null(n.obs)) {
      switch(fa, both = {
        if (sim) {
          legend("topright", c("  PC  Actual Data", "  PC  Simulated Data", 
                               " PC  Resampled Data", "  FA  Actual Data", 
                               "  FA  Simulated Data", " FA  Resampled Data"), 
                 col = c("black", "black", "black", "black", "black", 
                         "red"), pch = c(4, NA, NA, 2, NA, NA), 
                 text.col = "black", lty = c("solid", "dotted", 
                                              "dashed", "solid", "dotted", "dashed"), 
                 merge = TRUE)
        } else {
          legend("topright", c("  PC  Actual Data", " PC  Resampled Data", 
                               "  FA  Actual Data", " FA  Resampled Data"), 
                 col = c("black", "black", "black", "black"), pch = c(4, 
                                                                NA, 2, NA, NA), text.col = "black", lty = c("solid", 
                                                                                                             "dashed", "solid", "dashed"), merge = TRUE, 
                 bg = "black")
        }
      }, pc = {
        legend("topright", c("  PC  Actual Data", "  PC  Simulated Data", 
                             " PC  Resampled Data"), col = c("black", "black", 
                                                             "black", "black", "black", "black"), pch = c(4, NA, 
                                                                                                   NA, 2, NA, NA), text.col = "black", lty = c("solid", 
                                                                                                                                                "dotted", "dashed", "solid", "dotted", "dashed"), 
               merge = TRUE)
      }, fa = {
        legend("topright", c("  FA  Actual Data", "  FA  Simulated Data", 
                             " FA  Resampled Data"), col = c("black", "black", 
                                                             "black", "black", "black", "black"), pch = c(4, NA, 
                                                                                                   NA, 2, NA, NA), text.col = "black", lty = c("solid", 
                                                                                                                                                "dotted", "dashed", "solid", "dotted", "dashed"), 
               merge = TRUE)
      })
    }
    else {
      switch(fa, both = {
        legend("topright", c("PC  Actual Data", " PC  Simulated Data", 
                             "FA  Actual Data", " FA  Simulated Data"), 
               col = c("black", "black", "black", "black"), pch = c(4, 
                                                              NA, 2, NA), text.col = "black", lty = c("solid", 
                                                                                                       "dotted", "solid", "dotted"), merge = TRUE, 
               bg = "black")
      }, pc = {
        legend("topright", c("PC  Actual Data", " PC  Simulated Data"), 
               col = c("black", "black", "black", "black"), pch = c(4, 
                                                              NA, 2, NA), text.col = "black", lty = c("solid", 
                                                                                                       "dotted", "solid", "dotted"), merge = TRUE, 
               bg = "black")
      }, fa = {
        legend("topright", c("FA  Actual Data", " FA  Simulated Data"), 
               col = c("black", "black", "black", "black"), pch = c(4, 
                                                              NA, 2, NA), text.col = "black", lty = c("solid", 
                                                                                                       "dotted", "solid", "dotted"), merge = TRUE, 
               bg = "black")
      })
    }
  }
  abline(h = 1)
  if (fa != "pc") {
    abline(h = 0)
  }
  if (fa == "pc") {
    results <- list(fa.values = fa.valuesx, pc.values = valuesx, 
                    pc.sim = sim.pc, Call = cl)
    fa.test <- NA
  }
  else {
    results <- list(fa.values = fa.valuesx, fa.sim = sim.fa, 
                    pc.values = valuesx, pc.sim = sim.pc, Call = cl)
    if (sim) {
      fa.test <- which(!(fa.valuesx > sim.fa))[1] - 1
    }
    else {
      fa.test <- which(!(fa.valuesx > sim.far))[1] - 1
    }
    results$nfact <- fa.test
  }
  if (sim) {
    pc.test <- which(!(valuesx > sim.pc))[1] - 1
  }
  else {
    pc.test <- which(!(valuesx > sim.pcr))[1] - 1
  }
  results$ncomp <- pc.test
  cat("Parallel analysis suggests that ")
  cat("the number of factors = ", fa.test, " and the number of components = ", 
      pc.test, "\n")
  class(results) <- c("psych", "parallel")
  return(invisible(results))
}

#################################################################################################################
################################################ Modified Scree plot ############################################

scree1 <- function (rx, factors = TRUE, pc = TRUE, main = "Scree plot", 
                    hline = NULL, add = FALSE, use) 
{
  cl <- match.call()
  nvar <- dim(rx)[2]
  if (nvar != dim(rx)[1]) {
    rx <- cor(rx, use = use)
  }
  if (pc) {
    values <- eigen(rx)$values
    if (factors) {
      ylab = "Eigen values of factors and components"
      xlab = "factor or component number"
    }
    else {
      ylab = "Eigen values of components"
      xlab = " component number"
    }
  }
  else {
    values <- fa(rx)$values
    ylab = "Eigen values of factors"
    xlab = " factor number"
    factors <- FALSE
  }
  max.value <- max(values)
  if (!add) {
    plot(values, type = "b", main = main, pch = 16, ylim = c(0, 
                                                             max.value), ylab = ylab, xlab = xlab)
  }
  else {
    points(values, type = "b", , pch = 16)
  }
  if (factors) {
    fv <- fa(rx)$values
    points(fv, type = "b", pch = 21, lty = "dotted")
  }
  else {
    fv <- NULL
  }
  if (is.null(hline)) {
    abline(h = 1)
  }
  else {
    abline(h = hline)
  }
  if (factors & pc) {
    legend("topright", c("PC ", "FA"), pch = c(16, 21), text.col = "black", 
           lty = c("solid", "dotted"), merge = TRUE)
  }
  if (pc) {
    results <- list(fv = fv, pcv = values, call = cl)
  }
  else {
    results <- list(fv = values, pcv = NULL, call = cl)
  }
  class(results) <- c("psych", "scree")
  invisible(results)
}