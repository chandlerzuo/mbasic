#' @name MBASIC.MADBayes
#' @title MAD-Bayes method to fit the MBASIC model.
#' @param Y An N by I matrix containing the data from N experiments across I observation units (loci).
#' @param Gamma An N by I matrix for the prior estimated mean for the background state, for N experiments across the I observation units (loci).
#' @param fac A vector of length N denoting the experimental condition for each replicate.
#' @param lambdap,lambdaw,lambda Tuning parameters.
#' @param family The distribution of family to be used. Either "lognormal" or "negbin". See details for more information.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for error in checking the E-M algorithm's convergence. Default: 1e-04.
#' @param S The number of different states.
#' @param zeta The initial value of the proportion of unclustered units. Default: 0.2.
#' @param verbose Boolean variable for whether the model fitting messages are printed.
#' @param para A list of true paramters.
#' @details
#' TODO.
#' @useDynLib MBASIC
#' @return A list object.
#' @import cluster
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @export
MBASIC.MADBayes <- function(Y, Gamma, fac, lambdap = 0.5, lambdaw = 0.2, lambda = 5, maxitr = 100, S = 2, tol = 0.01, zeta = 0.1, verbose = TRUE, para = NULL) {
  ## Initialize
  ## prespecified
  K <- length(unique(fac))
  I <- ncol(Y)
  N <- nrow(Y)
  if(length(fac) != N)
    stop("Error: total number of replicates do not match with the number of rows in Y")

  if(prod(dim(Y) == dim(Gamma)) != 1)
    stop("Error: dimensions for Y and Gamma must be the same.")

  ## design matrix is N by K
  designMat <- matrix(0, ncol = K, nrow = N)
  for(k in 1:K) {
    designMat[fac == unique(fac)[k], k] <- 1
  }

  ## Scale the data from different replicates
  scaleFactor <- apply(Y, 1, mean)
  
  if(is.null(Gamma)) {
    Gamma <- Y - Y + 1
  }
  ## normalize the Gamma
  Gamma <- Gamma / rep(apply(Gamma, 1, mean), ncol(Gamma))
  Y <- Y / scaleFactor
  Gamma <- rbind(Gamma, matrix(0, ncol = ncol(Gamma), nrow = nrow(Gamma) * (S - 1)))
  for(s in seq(S)[-1]) {
    Gamma[(s - 1) * N + seq(N), ] <- rep(apply(Gamma[seq(N), ], 1, mean), I)
  }
  
  ## Initialize
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
  
  ## Initialize cluster
  if(verbose)
    message("Initialize clusters...")
  J <- max(c(as.integer(sqrt(I) / 4), 2))
  b <- sample(c(0, 1), I, prob = c(1 - zeta, zeta), replace = TRUE)
  clusterLabels <- sample(seq(J), I, replace = TRUE) - 1
  
  D <- apply(designMat, 1, function(x) which(x == 1)) - 1
  storage.mode(D) <- "integer"
  storage.mode(Theta) <- "integer"
  zeta <- mean(b)

  ret <- .Call("MADBayes", b, clusterLabels, Theta, Mu, D, Gamma, Y, lambdap, lambdaw, lambda, package = "MBASIC")

  allloss <- NULL
  allnclusters <- NULL

  associationMatrix <- matrix(0, nrow = I, ncol = I)
  maxId <- 0
  t0 <- Sys.time()
  mixed <- FALSE
  outliers <- rep(0, I)
  if(verbose)
    message("start iteration...")
  
  for(itr in seq(maxitr)) {
    ret <- .Call("MADBayes", ret$b, ret$clusterLabels, ret$Theta, ret$Mu, D, Gamma, Y, lambdap, lambdaw, lambda, package = "MBASIC")
    allloss <- c(allloss, ret$loss)
    allnclusters <- c(allnclusters, max(ret$clusterLabels) + 1)
    if(verbose & itr %% 10 == 0) {
      message(Sys.time() - t0, " have passed, number of iterations = ", itr)
    }
    if(length(unique(ret$clusterLabels)) != max(ret$clusterLabels) + 1)
      stop("Number of states is not equivalent to the maximum state index.")
    if(length(allloss) > 2)
      if(abs(diff(tail(allloss, 2))) < tol)
        break
  }

  if(maxitr <= itr)
    warning("MADBAYES procedure not converged.")

  for(j in unique(ret$clusterLabels)) {
    associationMatrix[ret$clusterLabels == j, ret$clusterLabels == j] <- 1
  }

  J <- max(ret$clusterLabels) + 1
  Z <- matrix(0, nrow = I, ncol = J)
  Z[cbind(seq(I), ret$clusterLabels + 1)] <- 1

  W <- ret$W[, seq(max(ret$clusterLabels) + 1)]
  
  Theta.err <- W.err <- ari <- mcr <- NULL
  if(!is.null(para)) {
    Theta.err <- sqrt(2 * sum(para$Theta != (ret$Theta + 1)) / I / K / S)
    W.f <- matrix(0, nrow = K * S, ncol = J)
    for(s in seq_len(S))
      W.f[ s + S * seq(0, K - 1), ] <- W[ seq_len(K) + K * (s - 1), ]
    mc <- matchCluster(W.f, para$W, Z, para$Z, ret$b, para$non.id)
    W.err <- mc$W.err
    ari <- mc$ari
    mcr <- mc$mcr
  }

  if(J == I) {
    sil.theta <- 0
  } else {
    sil.theta <- mean(silhouette(ret$clusterLabels, dist(t(ret$Theta), method = "manhattan"))[, 3])
  }

  ## Compute the loss of each term
  Theta.aug <- P.aug <- W.aug <- matrix(0, nrow = K * S, ncol = I)
  for(s in seq(S)) {
    Theta.aug[seq(K) + K * (s-1), ] <- as.integer(ret$Theta == s - 1)
    P.aug[seq(K) + K * (s-1), ] <- rep(ret$P[, s], each = K)
  }
  W.aug <- tcrossprod(W, Z)
  loss.p <- mean(abs(Theta.aug - P.aug)[, ret$b == 1])
  loss.w <- mean(abs(Theta.aug - W.aug)[, ret$b == 0])

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
  for(s in seq(S)) {
    denY <- dnorm(Y, mean = rep(Mu[, s], I), sd = rep(sqrt(Sigma[, s]), I), log = TRUE)
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
  dist.norm <- dist(t(Theta.norm), method = "manhattan")
  sil.norm <- mean(silhouette(ret$clusterLabels, dist.norm)[, 3])
  
  ## compute sihouette score from composite distance
  if(J == I) {
    sil.comp <- 0
  } else {
    dist.Y <- dist(t(Y - Mu.Y))
    dist.Theta <- dist(t(Theta.aug), method = "manhattan")
    sil.comp <- mean(silhouette(ret$clusterLabels, dmatrix = as.matrix(dist.Y) + lambdaw * as.matrix(dist.Theta))[, 3])
  }
  
  new("MBASICFit",
      Theta = t(ret$Theta) + 1,
      W = W,
      P = ret$P,
      b = ret$b,
      alllik = allloss,
      Mu = ret$Mu * scaleFactor,
      converged = (itr <= maxitr),
      Z = Z,
      AssociationMatrix = associationMatrix,
      Iter = itr,
      Theta.err = Theta.err,
      W.err = W.err,
      ARI = ari,
      MisClassRate = mcr,
      Loss = list(
        Y = loss.y,
        W = loss.w,
        P = loss.p,
        Silhouette.theta = sil.theta,
        Silhouette.norm = sil.norm,
        Silhouette.comp = sil.comp)
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
#' @param zeta The initial value of the proportion of unclustered units. Default: 0.2.
#' @param ncore The number of CPUs to be used for parallelization.
#' @param nfits The number of random restarts of the model.
#' @details
#' TODO.
#' @useDynLib MBASIC
#' @return A list object.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @import doMC
#' @export
MBASIC.MADBayes.full <- function(Y, Gamma, fac, lambdap = 15, lambdaw = 0.5, lambda = 20, maxitr = 100, S = 2, tol = 0.01, zeta = 0.1, ncore = 8, nfits = 1000, para = NULL) {
  require(doMC)
  registerDoMC(ncore)
  results <- foreach(i = seq(ncore)) %dopar% {
    set.seed(i + Sys.time())
    bestFit <-
      MBASIC.MADBayes(Y, Gamma, fac, lambdap = lambdap, lambdaw = lambdaw, lambda = lambda, maxitr = maxitr, S = S, tol = tol, zeta = zeta, verbose = FALSE, para = para)
    allLoss <- tail(bestFit@alllik, 1)
    allIter <- bestFit@Iter
    for(i in seq(as.integer(nfits / ncore))[-1]) {
      fit <-
        MBASIC.MADBayes(Y, Gamma, fac, lambdap = lambdap, lambdaw = lambdaw, lambda = lambda, maxitr = maxitr, S = S, tol = tol, zeta = zeta, verbose = FALSE, para = para)
      if(tail(fit@alllik, 1) < tail(bestFit@alllik, 1)) {
        bestFit <- fit
      }
      allLoss <- c(allLoss, tail(fit@alllik, 1))
      allIter <- c(allIter, bestFit@Iter)
    }
    list(BestFit = bestFit,
         Iter = allIter,
         Loss = allLoss)
  }
  bestFit <- results[[1]]$BestFit
  allIter <- results[[1]]$Iter
  allLoss <- results[[1]]$Loss
  for(i in seq(ncore)[-1]) {
    if(tail(bestFit@alllik, 1) > tail(results[[i]]$BestFit@alllik, 1)) {
      bestFit <- results[[i]]$BestFit
    }
    allIter <- c(allIter, results[[i]]$Iter)
    allLoss <- c(allLoss, results[[i]]$Loss)
  }
  return(list(BestFit = bestFit,
              Iter = allIter,
              Loss = allLoss)
       )
}
