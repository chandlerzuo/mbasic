#' @name MBASIC.MADBayes
#' @title MAD-Bayes method to fit the MBASIC model.
#' @param Y An N by I matrix containing the data from N experiments across I observation units (loci).
#' @param Mu0 An N by I matrix for the prior estimated mean for the background state, for N experiments across the I observation units (loci).
#' @param fac A vector of length N denoting the experimental condition for each replicate.
#' @param lambdap,lambdaw,lambda Tuning parameters.
#' @param family The distribution of family to be used. Either "lognormal" or "negbin". See details for more information.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for error in checking the E-M algorithm's convergence. Default: 1e-04.
#' @param S The number of different states.
#' @param zeta The initial value of the proportion of unclustered units. Default: 0.2.
#' @param verbose Boolean variable for whether the model fitting messages are printed.
#' @details
#' TODO.
#' @useDynLib MBASIC
#' @return A list object.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @export
MBASIC.MADBayes <- function(Y, Mu0, fac, lambdap = 0.5, lambdaw = 0.2, lambda = 5, maxitr = 100, S = 2, tol = 0.01, zeta = 0.1,
                            verbose = TRUE) {

  ## Initialize
  ## prespecified
  K <- length( unique( fac ) )
  I <- ncol( Y )
  N <- nrow( Y )
  if( length( fac ) != N )
    stop( "Error: total number of replicates do not match with the number of rows in Y" )

  if( prod( dim( Y ) == dim( Mu0 ) ) != 1 )
    stop( "Error: dimensions for Y and Mu0 must be the same." )

  ## design matrix D is K by N
  Dmat <- matrix( 0, nrow = K, ncol = length( fac ) )
  for( k in 1:K ){
    Dmat[ k, fac == unique( fac )[ k ] ] <- 1
  }
  ## normalize the Mu0
  Mu0 <- Mu0 * rep( apply( log( Y + 1 ), 1, mean ) / apply(Mu0, 1, mean ), ncol( Mu0 ) )
  Gamma <- t(Mu0)
  Y <- t(log(Y+1))
  Gamma <- cbind(Gamma, matrix(0, nrow = nrow(Gamma), ncol = ncol(Gamma) * (S - 1)))
  for(s in seq(S)[-1]) {
    Gamma[, (s - 1) * N + seq(N)] <- rep(apply(Gamma[, seq(N)], 2, mean), each = I)
  }
  ## Initialize
  Mu <- Sigma <- matrix(0, nrow = N, ncol = S)
  Theta <- matrix(0, nrow = I, ncol = K)
  storage.mode(Theta) <- "integer"
  foldChange <- Y / Gamma[, seq(N)]
  foldChange[Gamma[, seq(N)] == 0] <- 1
  avgFoldChange <- tcrossprod(foldChange, Dmat) / rep(apply(Dmat, 1, sum), each = I)
  for(k in seq(K)) {
    for(s in seq(S, 1)) {
      Theta[rank(avgFoldChange[, k]) <= s / S * I, k] <- s - 1
    }
  }
  DTheta <- tcrossprod(Theta, t(Dmat))
  for(n in seq(N)) {
    for(s in seq(S)) {
      Mu[n, s] <- mean(foldChange[DTheta[, n] == s - 1, n])
      Sigma[n, s] <- var(foldChange[DTheta[, n] == s - 1, n])
    }
  }
  ## Initialize cluster
  if(verbose)
    message("Initialize clusters...")
  J <- max(c(as.integer(sqrt(I) / 4), 2))
  if(FALSE) {
    d <- dist( Theta, method = "manhattan" )
    mind <- apply( as.matrix(d), 1, function( x ) min( x[x>0] ) )
    thr <- quantile( mind, 1 - zeta )
    b <- as.integer(mind > thr)
    hcfit <- hclust( d )
    States <- cutree( hcfit, k = J ) - 1
  } else {
    b <- sample(c(0, 1), I, prob = c(1 - zeta, zeta), replace = TRUE)
    States <- sample(seq(J), I, replace = TRUE) - 1
  }
  D <- apply(Dmat, 2, function(x) which(x == 1)) - 1
  storage.mode(D) <- "integer"
  storage.mode(Theta) <- "integer"
  zeta <- mean(b)

  ret <- .Call("MADBayes", b, States, Theta, Mu, D, Gamma, Y, lambdap, lambdaw, lambda, package = "MBASIC")

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
    ret <- .Call("MADBayes", ret$b, ret$States, ret$Theta, ret$Mu, D, Gamma, Y, lambdap, lambdaw, lambda, package = "MBASIC")
    allloss <- c(allloss, ret$loss)
    allnclusters <- c(allnclusters, max(ret$States) + 1)
    if(verbose & itr %% 10 == 0) {
      message(Sys.time() - t0, " have passed, number of iterations = ", itr)
    }
    if(length(unique(ret$States)) != max(ret$States) + 1)
      stop("Number of states is not equivalent to the maximum state index.")
    if(length(allloss) > 2)
      if(abs(diff(tail(allloss, 2))) < tol)
        break
  }

  if(maxitr <= itr)
    warning("MADBAYES procedure not converged.")

  for(j in unique(ret$States)) {
    associationMatrix[ret$States == j, ret$States == j] <- 1
  }

  J <- max(ret$States) + 1
  Z <- matrix(0, nrow = I, ncol = J)
  Z[cbind(seq(I), ret$States + 1)] <- 1

  allnclusters
  
  new("MBASICFit",
      Theta = t(ret$Theta),
      W = ret$W[, seq(max(ret$States) + 1)],
      P = ret$P,
      b = ret$b,
      alllik = allloss,
      Mu = ret$Mu,
      converged = (itr <= maxitr),
      Z = Z,
      AssociationMatrix = associationMatrix,
      Iter = itr
      )

}

#' @name MBASIC.MADBayes.full
#' @title MAD-Bayes method to fit the MBASIC model.
#' @param Y An N by I matrix containing the data from N experiments across I observation units (loci).
#' @param Mu0 An N by I matrix for the prior estimated mean for the background state, for N experiments across the I observation units (loci).
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
#' @export
MBASIC.MADBayes.full <- function(Y, Mu0, fac, lambdap = 15, lambdaw = 0.5, lambda = 20, maxitr = 100, S = 2, tol = 0.01, zeta = 0.1,
                            ncore = 8, nfits = 1000) {
  require(doMC)
  registerDoMC(ncore)
  results <- foreach(i = seq(ncore)) %dopar% {
    set.seed(i + Sys.time())
    bestFit <-
      MBASIC.MADBayes(Y, Mu0, fac, lambdap = lambdap, lambdaw = lambdaw, lambda = lambda, maxitr = maxitr, S = S, tol = tol, zeta = zeta, verbose = FALSE)
    allLoss <- tail(bestFit$Loss, 1)
    allIter <- bestFit$Iter
    for(i in seq(as.integer(nfits / ncore))[-1]) {
      fit <-
        MBASIC.MADBayes(Y, Mu0, fac, lambdap = lambdap, lambdaw = lambdaw, lambda = lambda, maxitr = maxitr, S = S, tol = tol, zeta = zeta, verbose = FALSE)
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
  return(results)
}
