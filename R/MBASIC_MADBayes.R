#' @name MBASIC.MADBayes
#' @title MAD-Bayes method to fit the MBASIC model.
#' @param Y An N by I matrix containing the data from N experiments across I observation units (loci).
#' @param Gamma An N by I matrix for the prior estimated mean for the background state, for N experiments across the I observation units (loci).
#' @param fac A vector of length N denoting the experimental condition for each replicate.
#' @param lambdaw,lambda Tuning parameters.
#' @param family The distribution of family to be used. Either 'lognormal' or 'negbin'. See details for more information.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for error in checking the E-M algorithm's convergence. Default: 1e-04.
#' @param S The number of different states.
#' @param verbose Boolean variable for whether the model fitting messages are printed.
#' @param para A list of true paramters.
#' @details
#' TODO.
#' @useDynLib MBASIC
#' @return A list object.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @export
MBASIC.MADBayes <- function(Y, Gamma, fac, lambdaw = 0.2, lambda = 200, maxitr = 100, S = 2, tol = 1e-6, verbose = TRUE, para = NULL, initialize = "kmeans") {

  ## Normalize Data
  NormalizeData()
  fit <- MBASIC.MADBayes.internal(Y = Y, Gamma = Gamma, fac = fac, lambdaw = lambdaw, lambda = lambda, maxitr = maxitr, S = S, tol = tol, verbose = verbose, para = para, initialize = initialize, Theta.init = NULL, Mu.init = NULL, Sigma.init = NULL, clusterLabels.init = NULL, scaleFactor = scaleFactor)
  return(fit)
}

#' @importFrom cluster silhouette
MBASIC.MADBayes.internal <- function(Y, Gamma, fac, lambdaw = NULL, lambda, maxitr = 20, S, tol = 1e-8, verbose = TRUE, para = NULL, initialize = "kmeans", Theta.init = NULL, Mu.init = NULL, Sigma.init = NULL, clusterLabels.init = NULL, scaleFactor, J = NULL) {

  GetModelStructure()

  facNames <- as.character(unique(fac))
  facMap <- seq(K)
  names(facMap) <- facNames
  fac <- as.character(fac)
  fac <- facMap[as.character(fac)]
  
  ## Initialize Mu, Sigma, Theta

  if(is.null(Mu.init) | is.null(Sigma.init) | is.null(Theta.init)) {
    InitializeTheta.MADBayes()
  } else {
    Theta <- Theta.init
    Mu <- Mu.init
    Sigma <- Sigma.init
  }

  if(is.null(lambdaw)) {
    lambdaw <- 2
  }
  
  lambdaw <- lambdaw * mean(na.omit(Sigma))
  lambda <- lambda * lambdaw
  
  ## Initialize cluster
  if(is.null(clusterLabels.init)) {
    InitializeClusters.MADBayes()
  } else {
    clusterLabels <- clusterLabels.init
  }
  
  b <- rep(0, I)
  
  ##  zeta <- mean(b)

  ret <- .Call("madbayes", clusterLabels, Theta, Mu, D, Gamma, Y, lambdaw, lambda, maxitr, tol, package = "MBASIC")

  t0 <- Sys.time()
  if(verbose)
    message("Finished iterations.")

  if(maxitr <= ret$Iter) {
    warning("MADBAYES procedure not converged.")
  }

  J <- max(ret$clusterLabels) + 1
  Z <- matrix(0, nrow = I, ncol = J)
  Z[cbind(seq(I), ret$clusterLabels + 1)] <- 1

  W <- ret$W[, seq(max(ret$clusterLabels) + 1), drop = FALSE]
  
  Theta.err <- W.err <- ari <- mcr <- numeric(0)
  if(!is.null(para)) {
    Theta.err <- sqrt(2 * sum(para$Theta != (ret$Theta + 1)) / I / K / S)
    W.f <- matrix(0, nrow = K * S, ncol = J)
    for(s in seq_len(S))
      W.f[ s + S * seq(0, K - 1), ] <- W[ seq_len(K) + K * (s - 1), ]
    mc <- matchCluster(W.f, para$W, Z, para$Z, ret$b, para$non.id)
    W.err <- mc$W.err
    mcr <- mc$mcr
    ## recompute ARI
    trueLabels <- apply(para$Z, 1, which.max)
    trueLabels[para$non.id] <- max(trueLabels) + seq(length(para$non.id))
    ari <- adjustedRandIndex(ret$clusterLabels, trueLabels)
  }

  ## Compute the loss of each term
  Theta.aug <- W.aug <- matrix(0, nrow = K * S, ncol = I)
  for(s in seq(S)) {
    Theta.aug[seq(K) + K * (s-1), ] <- as.integer(ret$Theta == s - 1)
  }
  W.aug <- tcrossprod(W, Z)
  loss.w <- mean((Theta.aug - W.aug) ^ 2)
  
  ## compute the loss of data fitting
  Theta.Y <- designMat %*% ret$Theta
  Mu.Y <- Sigma.Y <- Y - Y
  for(s in seq(S)) {
    id <- which(Theta.Y == s - 1)
    Mu.Y[id] <- (Gamma[seq(N) + N * (s - 1), ] * rep(Mu[, s], I))[id]
    Sigma.Y[id] <- rep(Sigma[, s], I)[id]
  }
  loss.y <- mean((Y - Mu.Y) ^ 2)

  ## Add normalized data
  ## Y: I * N, each column has mean 1
  Theta.norm <- matrix(0, ncol = I, nrow = K * S)
  PDF <- matrix(0, ncol = I, nrow = N * S)
  for(s in seq(S)) {
    denY <- dnorm(Y, mean = rep(Mu[, s], I), sd = rep(sqrt(Sigma[, s]), I), log = TRUE)
    PDF[seq(N) + N * (s - 1), ] <- denY
    Theta.norm[seq(K) + K * (s - 1), ] <- crossprod(designMat, denY)
  }
  ## for any i k, max_s Theta.norm[i, k, s] = 0
  Theta.norm <- Theta.norm - t(matrix(rep(apply(matrix(t(Theta.norm), nrow = I * K), 1, max), S), nrow = I))
  ## avoid getting Inf exponenant
  Theta.norm[Theta.norm > 5] <- 5
  Theta.norm <- exp(Theta.norm)
  Theta.total <- matrix(0, ncol = I, nrow = K)
  for(s in seq(S)) {
    Theta.total <- Theta.total + Theta.norm[(s - 1) * K + seq(K), ]
  }
  Theta.norm <- Theta.norm / t(matrix(rep(t(Theta.total), S), nrow = I))
  if(J == I | var(ret$clusterLabels) == 0) {
    sil.norm <- 0
  } else {
    dist.norm <- dist(t(Theta.norm), method = "manhattan")
    sil.norm <- mean(cluster::silhouette(ret$clusterLabels, dist.norm)[, 3])
  }

  P <- matrix(1 / S, ncol = S, nrow = I)
  V <- matrix(1, nrow = N, ncol = S)
  probz <- apply(Z, 2, mean)
  loglik <- .Call("loglik", W, P, V, 1e-10, probz, PDF, fac - 1, seq(S) - 1, package = "MBASIC")
  npars <- ncol(Z) - 1 + prod(dim(W)) * (S - 1) / S + N * S * 2

  Theta <- ret$Theta + 1
  rownames(Theta) <- facNames
  rownames(W) <- rep(facNames, S)
  Mu <- ret$Mu * scaleFactor
  Sigma <- ret$Sigma * scaleFactor * scaleFactor
  rownames(Mu) <- rownames(Sigma) <- rownames(V) <- facNames[fac]
 
  new("MBASICFit",
      Theta = Theta,
      W = W,
      clustProb = cbind(0, Z),
      alllik = ret$loss,
      Mu = Mu,
      Sigma = Sigma,
      converged = (ret$Iter <= maxitr),
      Z = Z,
      Iter = ret$Iter,
      Theta.err = Theta.err,
      W.err = W.err,
      ARI = ari,
      MisClassRate = mcr,
      Loss = list(
        lambdaw = lambdaw,
        lambda = lambda,
        Y = loss.y,
        W = loss.w,
        Silhouette = sil.norm,
        loglik = loglik,
        bic = -2 * loglik + log(N * I) * npars)
    )
}

#' @name MBASIC.MADBayes.full
#' @title MAD-Bayes method to fit the MBASIC model.
#' @param Y An N by I matrix containing the data from N experiments across I observation units (loci).
#' @param Gamma An N by I matrix for the prior estimated mean for the background state, for N experiments across the I observation units (loci).
#' @param fac A vector of length N denoting the experimental condition for each replicate.
#' @param lambdap,lambdaw,lambda Tuning parameters.
#' @param family The distribution of family to be used. Either "lognormal" or "negbin". See details for more information.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for error in checking the E-M algorithm's convergence. Default: 1e-04.
#' @param S The number of different states.
#' @param ncores The number of CPUs to be used for parallelization.
#' @param nfits The number of random restarts of the model.
#' @details
#' TODO.
#' @useDynLib MBASIC
#' @return A list object including the following fields:
#' \tabular{ll}{
#' allFits \tab A list of \linkS4class{MBASICFit} objects for the best model fit with each lambda.\cr
#' lambda \tab A vector of all lambdas corresponding to \code{allFits}.\cr
#' Loss \tab A vector for the loss corresponding to \code{allFits}.\cr
#' BestFit \tab The \linkS4class{MBASICFit} object with largest Silhouette score.\cr
#' Iter \tab Number of iterations for \code{BestFit}.\cr
#' Time \tab Time in seconds used to fit the model.\cr
#' }
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @import foreach
#' @export
MBASIC.MADBayes.full <- function(Y, Gamma = NULL, fac, lambdaw = NULL, lambda = NULL, maxitr = 30, S = 2, tol = 1e-10, ncores = 15, nfits = 1, nlambdas = 30, para = NULL, initialize = "kmeans") {
  t0 <- Sys.time()
  if(!is.null(lambda)) {
    ncores <- min(c(ncores, length(lambda) * nfits))
  }
  
  startParallel(ncores)
  
  if(!is.null(lambda)) {
    lambdas <- sort(unique(lambda))
    alllambdas <- rep(lambdas, each = nfits)
    results <- foreach(i = seq_along(alllambdas)) %dopar% {
      set.seed(i + Sys.time())
      fit <- MBASIC.MADBayes(Y, Gamma, fac, lambdaw = lambdaw, lambda = alllambdas[i], maxitr = maxitr, S = S, tol = tol, verbose = FALSE, para = para, initialize = initialize)
      list(fit = fit, lambda = alllambdas[i])
    }
    initLosses <- NULL
  } else {
    if(!is.numeric(nlambdas)) {
      stop("Error: 'nlambdas' must take a numeric value.")
    }
    NormalizeData()
    GetModelStructure()
    ## Initialize Theta, Mu, Sigma
    InitializeTheta.MADBayes()
    ## Initialize clusters and pick a range of lambda values
    if(initialize == "madbayes") {
      ret <- foreach(i = seq(as.integer(sqrt(I)) + 1)) %dopar% {
        .Call("madbayes_init", Theta, 0, S, i, package = "MBASIC")
      }
    } else if(initialize == "kmeans++") {
      ret <- foreach(i = seq(as.integer(sqrt(I)) + 1)) %dopar% {
        .Call("madbayes_init_kmeanspp", Theta, S, i, package = "MBASIC")
      }
    } else if(initialize == "kmeans") {
      Theta.aug <- matrix(0, nrow = K * S, ncol = I)
      for(s in seq(S)) {
        Theta.aug[seq(K) + (s - 1) * K, ] <- (Theta == s)
      }
      ret <- foreach(i = seq(as.integer(sqrt(I)) + 1)) %dopar% {
        fit.kmeans <- kmeans(t(Theta.aug), i)
        return(list(loss = fit.kmeans$tot.withinss * 2, clusterLabels = fit.kmeans$cluster - 1))
      }
    } else {
      stop("Error: a vector for 'lambda' values must be provided.")
    }

    endParallel()
    
    message("Initialized clusters")
    
    initLosses <- rep(0, length(ret))
    allClusterLabels <- list()
    for(i in seq_along(ret)) {
      initLosses[i] = ret[[i]]$loss
      allClusterLabels[[i]] <- ret[[i]]$clusterLabels
    }

    slopes <- abs(diff(sort(initLosses, decreasing = TRUE)))
    slopes <- slopes[-1]
    allLambdas <- (slopes[-1] + slopes[-length(slopes)]) / 2
    allLambdas <- sort(allLambdas, decreasing = TRUE)
    minLambda <- min(allLambdas)
    maxLambda <- max(allLambdas)
    ## lambdas <- seq(minLambda, maxLambda, length = nlambdas)
    lambdas <- unique(quantile(allLambdas, seq(nlambdas) / (nlambdas + 2)))
    lambdas <- lambdas[-c(1, length(lambdas))]
    alllambdas <- rep(lambdas, each = nfits)

    initClusterLabels <- list()
    usedids <- numeric(0)
    freeids <- seq_along(allLambdas)
    allJs <- seq_along(alllambdas)
    for(i in seq_along(alllambdas)) {
      j <- freeids[which.min(abs(allLambdas[freeids] - alllambdas[i]))[1]]
      initClusterLabels[[i]] <- allClusterLabels[[j + 2]]
      usedids <- c(usedids, j)
      freeids <- setdiff(freeids, usedids)
      j <- which.min(abs(allLambdas - alllambdas[i]))[1]
      allJs[i] <- max(allClusterLabels[[j + 2]]) + 1
    }

    message("Selected lambda values: ", paste(allLambdas, collapse = ", "))
    
    results <- foreach(i = seq_along(alllambdas)) %dopar% {
      set.seed(i + Sys.time())
      ##      fit <- MBASIC.MADBayes.internal(Y, Gamma, fac, lambdaw = lambdaw, lambda = alllambdas[i], maxitr = maxitr, S = S, tol = tol, verbose = FALSE, para = para, initialize = initialize, Theta.init = Theta, Mu.init = Mu, Sigma.init = Sigma, clusterLabels.init = initClusterLabels[[i]], scaleFactor = scaleFactor)
      fit <- MBASIC.MADBayes.internal(Y, Gamma, fac, lambdaw = lambdaw, lambda = alllambdas[i], maxitr = maxitr, S = S, tol = tol, verbose = TRUE, para = para, initialize = initialize, Theta.init = Theta, Mu.init = Mu, Sigma.init = Sigma, scaleFactor = scaleFactor, J = allJs[i])
      list(fit = fit, lambda = alllambdas[i])
    }
  }

  if(.Platform$OS.type != "unix") {
    stopCluster(cl)
  }
  
  message("Finished individual models")
  ## Within the same lambda, choose the model that minimizes the loss
  bestLosses <- rep(Inf, length(lambdas))
  bestFits <- list()
  bestIters <- rep(0, length(lambdas))
  for(i in seq_along(results)) {
    lambdaid <- which(lambdas == alllambdas[i])[1]
    if(bestLosses[lambdaid] > tail(results[[i]]$fit@alllik, 1)) {
      bestLosses[lambdaid] <- tail(results[[i]]$fit@alllik, 1)
      bestFits[[lambdaid]] <- results[[i]]$fit
      bestIters[lambdaid] <- results[[i]]$fit@Iter
    }
  }
  
  ## Between different lambdas, choose the model with the largest Silhouette score
  bestSil <- -Inf
  bestFit <- NULL
  for(fit in bestFits) {
    if(fit@Loss$Silhouette > bestSil) {
      bestFit <- fit
      bestSil <- fit@Loss$Silhouette
    }
  }
  bestbic <- Inf
  bestFit.bic <- NULL
  for(fit in bestFits) {
    if(fit@Loss$bic < bestbic) {
      bestFit.bic <- fit
      bestbic <- fit@Loss$bic
    }
  }
  
  return(list(allFits = bestFits,
              BestFit = bestFit,
              BestFit.bic = bestFit.bic,
              Iter = bestIters,
              Loss = bestLosses,
              lambda = lambdas,
              Time = as.numeric(Sys.time() - t0, units = "secs"),
              InitLoss = initLosses)
         )
}

InitializeTheta.MADBayes <- function() {
  Inherit(c("N", "S", "K", "I", "Y", "Gamma", "designMat", "D", "maxitr", "tol"))
  Mu <- Sigma <- matrix(0, nrow = N, ncol = S)
  Theta <- matrix(0, nrow = K, ncol = I)
  storage.mode(Theta) <- "integer"
  foldChange <- Y / Gamma[seq(N), ]
  foldChange[Gamma[seq(N), ] == 0] <- 1
  avgFoldChange <- crossprod(foldChange, designMat) / rep(apply(designMat, 2, sum), each = I)
  for(k in seq(K)) {
    for(s in seq(S, 1)) {
      Theta[k, rank(avgFoldChange[, k]) <= s / S * I] <- s - 1
    }
  }
  ## DTheta: N by I
  DTheta <- designMat %*% Theta
  for(n in seq(N)) {
    for(s in seq(S)) {
      Mu[n, s] <- mean(foldChange[n, DTheta[n, ] == s - 1])
      Sigma[n, s] <- var(foldChange[n, DTheta[n, ] == s - 1])
    }
  }
  Sigma[Sigma <= 0] <- min(Sigma[Sigma > 0])
  
  storage.mode(Theta) <- "integer"
  
  ## initialize Theta
  
  ret <- .Call("madbayes_theta", Theta, Mu, D, Gamma, Y, maxitr, tol, package = "MBASIC")

  Theta <- ret$Theta
  Mu <- ret$Mu
  Sigma <- ret$Sigma

  Return(c("Theta", "Mu", "Sigma"))

}

InitializeClusters.MADBayes <- function() {
  Inherit(c("verbose", "I", "initialize", "K", "S", "Theta", "lambda", "lambdaw", "J"))
  if(verbose)
    message("Initialize clusters...")
  
  ## Sample J from an exponential distribution with the median sqrt(I)/4
  if(is.null(J)) {
    J <- sample(seq(I)[-1], 1, prob = exp(-seq(I)/sqrt(I)*4*log(2))[-1])
  } else  {
    J <- as.integer(J)
  }

  if(initialize == "kmeans") {
    ## use K-means to initialize clusters
    Theta.aug <- matrix(0, nrow = K * S, ncol = I)
    for(s in seq(S)) {
      Theta.aug[seq(K) + (s - 1) * K, ] <- (Theta == s)
    }
    ##  clusterLabels <- sample(seq(J), I, replace = TRUE) - 1
    clusterLabels <- kmeans(t(Theta.aug), centers = J)$cluster - 1
  } else if(initialize == "kmeans++" ) {
    ret <- .Call("madbayes_init_kmeanspp", Theta, S, J, package = "MBASIC")
    clusterLabels <- ret$clusterLabels
  } else if(initialize == "madbayes") {
    ret <- .Call("madbayes_init", Theta, lambda / lambdaw, S, I, package = "MBASIC")
    clusterLabels <- ret$clusterLabels
  } else {
    clusterLabels <- sample(seq(J), I, replace = TRUE)
  }
  Return("clusterLabels")
}

GetModelStructure <- function() {

  Inherit(c("fac", "Y", "Gamma"))
  ## prespecified
  K <- length(unique(fac))
  I <- ncol(Y)
  N <- nrow(Y)
  if(length(fac) != N)
    stop("Error: total number of replicates do not match with the number of rows in Y")

  ## design matrix is N by K
  designMat <- matrix(0, ncol = K, nrow = N)
  for(k in 1:K) {
    designMat[fac == unique(fac)[k], k] <- 1
  }
  D <- apply(designMat, 1, function(x) which(x == 1)) - 1
  storage.mode(D) <- "integer"

  Return(c("K", "I", "N", "designMat", "D"))
}

NormalizeData <- function() {
  Inherit(c("Y", "Gamma", "S"))
  ## Scale the data from different replicates
  scaleFactor <- apply(Y, 1, mean)
  N <- length(scaleFactor)
  I <- ncol(Y)
  
  if(is.null(Gamma)) {
    Gamma <- matrix(1, nrow = nrow(Y), ncol = ncol(Y))
  }
  
  if(prod(dim(Y) == dim(Gamma)) != 1) {
    stop("Error: dimensions for Y and Gamma must be the same.")
  }
  
  ## normalize the Gamma
  Gamma <- Gamma / rep(apply(Gamma, 1, mean), ncol(Gamma))
  Y <- Y / scaleFactor
  Gamma <- rbind(Gamma, matrix(0, ncol = ncol(Gamma), nrow = nrow(Gamma) * (S - 1)))
  for(s in seq(S)[-1]) {
    Gamma[(s - 1) * N + seq(N), ] <- rep(apply(Gamma[seq(N), ], 1, mean), I)
  }

  Return(c("Y", "Gamma", "scaleFactor"))

}
