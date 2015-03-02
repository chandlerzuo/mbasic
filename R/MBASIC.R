#' @name MBASIC
#' @title Bayesian clustering model for a state-space matrix.
#'
#' @description This function is designed to analyze general state-space models. The data consists of observations over I units under N experiments with K different conditions. There are S states for each experiment and unit.
#' @param Y An N by I matrix containing the data from N experiments across I observation units.
#' @param S An integer for the number of states.
#' @param fac A vector of levels repr1esenting the conditions of each replicate.
#' @param struct A K by J matrix indicating the structures of each cluster.
#' @param J The number of clusters to be identified.
#' @param family The distribution of family to be used. Either "lognormal", "negbin", "binom", "gamma-binom". See details for more information.
#' @param method A string for the fitting method, 'MBASIC' (default), 'PE-MC', 'SE-HC',  or 'SE-MC'. See details for more information.
#' @param para A list object that contains the true model parameters. Default: NULL. See details for more information.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for error in checking the E-M algorithm's convergence. Default: 1e-04.
#' @param zeta The initial value for the proportion of units that are not clustered. Default: 0.1. If 0, no singleton cluster is fitted.
#' @param out The file directory for writing fitting information in each E-M iteration. Default: NULL (no information is outputted).
#' @param verbose A boolean variable indicating whether intermediate model fitting metrics should be printed. Default: FALSE.
#' @param statemap A vector the same length as the number of mixture components, and taking values from 1 to S representing the states of each component. Default: NULL.
#' @details
#' Function MBASIC currently supports two different distributional families: log-normal and negative binomial. This should be specified by the 'family' argument.\cr
#' For the log-normal distributions, log(Y+1) is modeled as normal distributions. For experiment n, if locus i has state s, distribution for log(Y[n,i]+1) is N(Mu[n,s], Sigma[n,s]).\cr
#' For the negative binomial distributions, the meanings of Mu and Sigma are different. For experiment n, if locus i has state s, distribution of Y[n,i] is NB(Mu[n,s], Sigma[n,s]). In this package, NB(mu, a) denotes the negative-binomial distribution with mean mu and size a (i.e. the variance is mu*(1+mu/a)).\cr
#'  The 'method' argument determines what fitting method will be used. The default is 'MBASIC', where the states and the clustering are simultaneously estimated. 'SE-HC' and 'SE-MC' methods use 2-step algorithms. In the first step, both estimate the states for each unit by an E-M algorithm for each experiment. In the second step, 'SE-HC' uses hierarchical clustering to cluster the units, while 'SE-MC' uses function 'MBASIC.state' to identify clusters.\cr
#' The 'para' argument takes a list object that is supposed to include the following fields:
#'\tabular{ll}{
#' W \tab A K by (J*S) matrix. The (k,J*(s-1)+j)-th entry is the probability that the units in cluster j has state s in the k-th experiment.\cr
#' Z \tab An I by J matrix. The (i,j)-th entry is the indicator whether the i-th unit belongs to cluster j.\cr
#' Theta \tab A K by (I*S) matrix. The (k,I*(s-1)+i)-th entry is the probability that the i-th unit has state s in the k-th experiment.\cr
#' non.id \tab A binary vector of length I. The i-th entry is the indicator whether the i-th unit does not belong to any cluster.
#' }
#' This argument is intended to carry the true parameters in simulation studies. If it is not null, then the model also computes a number of metrics that describes the error in model fitting. Users should be cautious that the order of the rows and columns of matrices in the fields of para should match the Y matrix.
#' @return An object of class 'MBASICFit'.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' ## Simulate a dataset
#' dat.sim <- MBASIC.sim(xi = 2, family = "lognormal", I = 1000, fac = rep(1:10, each = 2), J = 3, S = 3, zeta = 0.1)
#' ## Fit the model
#' dat.sim.fit <- MBASIC(Y = dat.sim$Y, S = 3, fac = rep(1:10, each = 2), J = 3, maxitr = 3, para = NULL, family = "lognormal", method = "MBASIC", zeta = 0.1, tol = 1e-04)
#' @useDynLib MBASIC
#' @export
MBASIC <- function(Y, S, fac, J=NULL, maxitr = 100, struct = NULL, para = NULL,  family="lognormal", method = "MBASIC", zeta = 0.1, tol = 1e-4, out = NULL, X = NULL, verbose = FALSE, statemap = NULL) {

  write.out(out, "Started")
  if(! method %in% c("SE-HC", "SE-MC", "PE-MC", "MBASIC")) {
    message("Error: 'method' must be one of 'SE-HC', 'SE-MC', 'PE-MC' or 'MBASIC'.")
    return
  }

  if(! family %in% c("lognormal", "negbin", "binom", "gamma-binom", "scaled-t")) {
      message("Error: 'family' must be one of 'lognormal', 'negbin', 'binom', 'gamma-binom', 'scaled-t'.")
  }

  if(is.null(statemap)) {
    ## states and components have one-to-one matching
    statemap <- seq(S)
  }

  M <- length(statemap)
  stateMap <- matrix(0, nrow = M, ncol = S)
  stateMap[cbind(seq(M), statemap)] <- 1
  
  ## prespecified
  K <- length(unique(fac))
  I <- ncol(Y)
  N <- nrow(Y)
  if(length(fac) != N)
    message("Error: total number of replicates do not match with the number of rows in Y")
  
  if(is.null(struct)) {
    if(is.null(J))
      message("Error: either struct or J must not be missing.")
    struct <- matrix(seq_len(K), nrow = K, ncol = J)
  } else {
    if(is.null(J))
      J <- ncol(struct)
    J <- sum(J)
    if(ncol(struct)!= sum(J) | nrow(struct) != K)
      message("Error: the dimension of struct is inconsistent with grouping structure!")
  }

  if(prod(sort(unique(statemap)) == seq(S)) != 1) {
    stop("Error: 'statemap' must consists of values 1 to 'S'.")
  }

  ## design matrix D is K by N
  designMap <- matrix(0, nrow = N, ncol = K)
  for(k in seq(K)) {
    designMap[fac == unique(fac)[k], k] <- 1
  }
  SampleToExp <- apply(designMap, 1, function(x) which(x==1))
  unitMap <- matrix(0, nrow = N * M, ncol = N * S)
  for(m in seq(M)) {
    rows <- (m - 1) * N + seq(N)
    cols <- (statemap[m] - 1) * N + seq(N)
    unitMap[cbind(rows, cols)] <- 1
  }

  numpar <- I * (S - 1) + ## P
    J + ## probz and zeta
      sum(apply(struct, 2, function(x) length(unique(x)))) * (S - 1) + ## W
        N * M + ## Mu
          N * (M - S) ## V

  if(family != "binom") {
    numpar <- numpar + N * M ## Sigma
  }
  
  outitr <- 0
  totallik <- oldlik <- 0
  alllik <- allerr <- allzeta <- bestW <- bestV <- allmisclass <- matchId1 <- W.err <- matchId2 <- allari <- numeric(0)
  maxlik <- -Inf
  
  write.out(out, "Initialized parameters")
  
  ## initialize distributions
  V <- Sigma <- Mu <- matrix(0, nrow = N, ncol = M)

  InitDist()
  
  b <- rep(0, I)
  
  ## initialize the matrices by hierarchical clustering
  ## in constructing Z, cluster all locis
  ## This gives deterministic initialization
  ProbMat <- matrix(0, nrow = K * S, ncol = I)
  InitStates()

  if(!is.null(para)) {
    ProbMat.true <- matrix(0, nrow = S * K, ncol = I)
    for(s in seq(S)) {
        ProbMat.true[(s - 1) * K + seq(K), ] <- as.numeric(para$Theta == s)
    }
  }

  if(method != "MBASIC") {
    ## SE-HC, PE-MC or SE-MC method
    
    ## em step to estimate Theta
    allpar <- c(c(V), c(Mu), c(Sigma), c(Pi))
    oldpar <- 0
    for(itr in 1:maxitr) {
      ## check for convergence
      if(max(abs(oldpar - allpar))< tol)
        break
      
      ## M step
      UpdateDist()
      
      ## E step
      UpdateStates()
      oldpar <- allpar
      allpar <- c(c(V), c(Mu), c(Sigma), c(Pi))
    }## finish iteration
    
    Theta <- matrix(-1, nrow = K, ncol = I)
    for(k in seq_len(K)) {
      idx <- k + K * (seq_len(S) - 1)
      Theta[k,] <- apply(ProbMat[idx,], 2, which.max)
    }

    if(!is.null(para))
      allerr <- sqrt(sum((ProbMat - ProbMat.true) ^ 2) / I / K / (S - 1))
    
    if(method != "PE-MC") {
      ret <- MBASIC.state(Theta, J=J, zeta = zeta, struct = struct, method = method, maxitr = maxitr, tol = tol, para = para, out = out)
      
      conv <- FALSE
      if(ret@converged & itr < maxitr)
        conv <- TRUE
      
      ## Pi is the proportion for components in the k experiment to have state s
      ## Pi is different from Z. Z is the posterior probability.

      Mu.err <- numeric(0)
      if(prod(dim(Mu) == dim(para$Mu)) == 1) {
        Mu.err <- sqrt(mean((Mu - para$Mu) ^ 2))
      }

      write.out(out, paste("mis-class rate ", ret@MisClassRate))
      write.out(out, paste("Error for W ",  round(ret@W.err, 3)))
      write.out(out, paste("ARI ", ret@ARI))
      write.out(out, paste("loglik", round(tail(ret@alllik, 1), 3), "err", round(allerr, 3)))
      write.out(out, paste("Error for Mu", round(Mu.err, 3)))

      return(new("MBASICFit",
                 Theta = ProbMat,
                 W = ret@W,
                 Z = ret@Z,
                 V = V,
                 b = ret@b,
                 clustProb = ret@clustProb,
                 lik = ret@lik,
                 alllik = ret@alllik,
                 zeta = ret@zeta,
                 Mu = Mu,
                 Sigma = Sigma,
                 probz = ret@probz,
                 P = ret@P,
                 converged = conv,
                 Theta.err = allerr,
                 ARI = ret@ARI,
                 W.err = ret@W.err,
                 MisClassRate = ret@MisClassRate,
                 Mu.err = Mu.err
                 )
             )
    }
  }
  
  ## initialize W, Z, b
  ## ProbMat <- D.rep %*% ProbMat
  d <- dist(t(ProbMat))
  mind <- apply(as.matrix(d), 1, function(x) min(x[x>0]))
  thr <- quantile(mind, 1 - zeta)
  id <- which(mind < thr)
  b <- rep(1, I)
  b[id] <- 0
  d <- dist(t(ProbMat[,id]))
  fit <- hclust(d)
  groups <- cutree(fit, k = J)
  Z <- matrix(0, nrow = I, ncol = J)
  Z[cbind(1:I, sample(1:J, I, replace = TRUE))] <- 1
  Z[id,] <- 0
  Z[cbind(id, groups)] <- 1
  W <- matrix(1/S, nrow = K * S, ncol = J)
  for(j in seq_len(J)) {
    if(length(id[groups == j]) > 0) {
      W[, j] <- apply(t(ProbMat[ , id[groups == j]]), 2, mean)
    }
  }
  predZ <- Zcond <- Z
  b.prob <- b
  clustOrder <- .orderCluster(W, struct)
  W <- W[, clustOrder]
  Z <- Z[, clustOrder]
  W <- .structure(W, struct)
  ## initialize p, probz
  P <- matrix(0, nrow = I, ncol = S)
  for(s in seq_len(S)) {
    idx <- seq_len(K) + K * (s - 1)
    P[, s] <- apply(ProbMat[idx, ], 2, mean)
  }
  probz <- apply(rbind(Z, diag(rep(1, J))), 2, mean)
  
  ## EM algorithm
  ## change storage modes for C modules
  storage.mode(statemap) <- "integer"
  storage.mode(fac) <- "integer"
  oldpar <- 0
  newpar <- c(c(W), probz, zeta, c(P), c(V), c(Mu), c(Sigma))
  
  for(outitr in seq_len(maxitr)) {
    
    if(max(abs(newpar - oldpar)) < tol) {
      break
    }
    
    if(outitr == 1 | method == "MBASIC") {
        ## only compute the PDF once if the method is PE-MC
        PDF <- matrix(0, nrow = N * M, ncol = I)
        for(m in seq_len(M)) {
            PDF[seq(N) + (m - 1) * N, ] <- logdensity(Y, Mu[, m], Sigma[, m], X, family)
        }
        PDF <- trimLogValue(PDF)
    }
    
    oldlik <- totallik
    totallik <- .Call("loglik", W, P, V, zeta, probz, PDF, fac - 1, statemap - 1, package="MBASIC")
    if(verbose) {
        write.out(out, paste("itr", outitr, "lik", round(tail(totallik, 1), 2), "zeta", round(zeta, 2)))
    }
    
    alllik <- c(alllik, totallik)
    allzeta <- c(allzeta, zeta)
    Theta <- matrix(-1, nrow = K, ncol = I)
    for(k in seq_len(K)) {
        idx <- k + K * (seq_len(S) - 1)
        Theta[k,] <- apply(ProbMat[idx,], 2, which.max)
    }

    if(length(para) > 1 & verbose) {
        PrintUpdate()
    }
    
    if(maxlik < totallik) {
        maxlik <- totallik
        bestb <- b.prob
        bestW <- W
        bestTheta <- Theta
        bestV <- V
    }
    
    ## E step
    ## M step for some parameters
    estep.result <- .Call("e_step", W, P, V, zeta, probz, PDF, fac - 1, statemap - 1, package = "MBASIC")
    if(FALSE) {
        estep.r1 <- estep(W, P, V, zeta, probz, PDF, designMap, stateMap, unitMap)
        for(s in names(estep.result)) {
            print(s)
            print(max(abs(estep.result[[s]] - estep.r1[[s]])))
        }
    }
    
    ## Expected Theta matrix
    ProbMat <- estep.result[["Theta"]]
    ProbMat.full <- estep.result[["Theta_nu"]]
    ## Maximizers
    zeta <- estep.result[["zeta"]]
    W <- estep.result[["W"]]
    probz <- estep.result[["probz"]]
    predZ <- estep.result[["Z"]]
    Zcond <- estep.result[["Zcond"]]
    b.prob <- estep.result[["b_prob"]]
    V <- estep.result[["V"]]
    
    clustOrder <- .orderCluster(W, struct)
    W <- W[, clustOrder]
    W <- .structure(W, struct)
    probz <- probz[clustOrder]
    Zcond <- Zcond[, clustOrder]
    predZ <- predZ[, clustOrder]
    
    if(method != "PE-MC") {
        ## skip this step if the method is PE-MC
        ## M-step for Mu and Sigma
        UpdateDist()
    }
    
    oldpar <- newpar
    newpar <- c(c(W), probz, zeta, c(P), c(V), c(Mu), c(Sigma))
    
  }## finish outer loop
  
  conv <- FALSE
  if(outitr < maxitr)
      conv <- TRUE
  
  if(length(para) > 1)
      PrintUpdate()
  
  new("MBASICFit",
      Theta = ProbMat,
      W = bestW,
      V = bestV     ,
      Z = predZ,
      b = bestb,
      clustProb = cbind(b, Zcond * b),
      aic = - 2 * tail(alllik, 1) + 2 * numpar,
      bic = - 2 * tail(alllik, 1) + log(N * I) * numpar,
      aicc = -2 * tail(alllik, 1) + 2 * numpar + 2 * numpar * (numpar + 1) / (N * I - numpar - 1),
      alllik = alllik,
      lik = tail(alllik, 1),
      zeta = zeta,
      Mu = Mu,
      Sigma = Sigma,
      probz = probz,
      P = P,
      converged = conv,
      Theta.err = allerr,
      ARI = tail(allari, 1),
      W.err = tail(W.err, 1),
      MisClassRate = tail(allmisclass, 1),
      Struct = struct,
      Mu.err = Mu.err
    )
}

InitStates <- function() {
    Inherit()
    ## initialize V
    for(s in seq(S)) {
      ids <- which(statemap == s)
      V[, ids] <- 1 / length(ids)
    }
    ## initialize ProbMat
    totalF <- matrix(0, nrow = K, ncol = I)
    F1 <- matrix(0, nrow = K * S, ncol = I)
    totalF.full <- matrix(0, nrow = N, ncol = I)
    F1.full <- matrix(0, nrow = N * M, ncol = I)
    for(m in seq(M)) {
      idx <- (m - 1) * N + seq_len(N)
      F1.full[idx,] <- logdensity(Y, Mu[, m], Sigma[, m], X, family)
    }
    F1.full <- trimLogValue(F1.full)
    F1.full <- exp(F1.full)
    ## convert into (NS) x I
    F1.tmp <- log(crossprod(unitMap, F1.full))
    
    ## initialize ProbMat.full
    for(m in seq(M)) {
      idx <- (m - 1) * N + seq(N)
      F1.full[idx, ] <- exp(F1.full[idx, ])
      totalF.full <- totalF.full + F1.full[idx, ]
    }
    totalF.full <- t(matrix(rep(c(t(totalF.full)), M), nrow = I))
    ProbMat.full <- F1.full / totalF.full
    ProbMat.full[totalF.full == 0] <- 1 / M
    ProbMat.full <- trimProbValue(ProbMat.full)
    
    ## convert to KS by I
    for(s in seq(S)) {
      idx <- (s - 1) * K + seq_len(K)
      F1[idx, ] <- exp(crossprod(designMap, F1.tmp[(s - 1) * N + seq_len(N), ]))
      totalF <- totalF + F1[idx, ]
    }
    totalF <- t(matrix(rep(c(t(totalF)), S), nrow = I))
    ProbMat <- F1 / totalF
    ProbMat[totalF == 0] <- 1/S
    ProbMat <- trimProbValue(ProbMat)

    Pi <- matrix(rep(apply(stateMap, 2, sum) / M, each = K), nrow = K, ncol = S)

    assign("ProbMat.full", ProbMat.full, envir = parent.frame())
    assign("ProbMat", ProbMat, envir = parent.frame())
    assign("V", V, envir = parent.frame())
    assign("Pi", Pi, envir = parent.frame())
}

InitDist <- function() {
    Inherit()
    for(m in seq(M)) {
        if(family == "lognormal") {
            Y.sec <- c(Y)[c(Y) <= quantile(c(Y), m / M) & c(Y) >= quantile(c(Y), (m - 1) / M)]
            Y.sec <- log(Y.sec + 1)
	} else if(family  == "negbin") {
            Y.sec <- c(Y)[c(Y) < quantile(c(Y), m / M) & c(Y) > quantile(c(Y), (m - 1) / M)]
        } else if(family == "scaled-t") {
	    Y.sec <- abs(Y)
            Y.sec <- c(Y.sec)[c(Y.sec) < quantile(c(Y.sec), m / M) & c(Y.sec) > quantile(c(Y.sec), (m - 1) / M)]
	} else if(family == "gamma-binom") {
            ## gamma-binomial distribution
            ratio <- Y / X
            ratio[X == 0] <- mean(Y[X > 0] / X[X > 0])
            Y.sec <- c(ratio)[ratio <= quantile(ratio, m / M) & ratio >= quantile(ratio, (m - 1) / M)]
        } else {
	    ratio <- Y / X
            id1 <- which(X > 0)
            ratio <- ratio[id1]
	    id <- which(ratio <= quantile(ratio, m / M) & ratio >= quantile(ratio, (m - 1) / M))
	    Y.sec <- sum(Y[id1[id]]) / sum(X[id1[id]])
	}
        m1 <- mean(Y.sec)
        m2 <- mean(Y.sec * Y.sec)
        MomentEstimate()
    }
    assign("Mu", Mu, envir = parent.frame())
    assign("Sigma", Sigma, envir = parent.frame())
}

UpdateDist <- function() {
    Inherit()
    
    for(m in seq_len(M)) {
        idx <- seq(N) + (m - 1) * N
        if(family == "lognormal") {
            m1 <- apply(log(Y + 1) * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
            m2 <- apply(log(Y + 1) * log(Y + 1) * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
	} else if(family == "negbin"){
            ## negative binomial family
            m1 <- apply(Y * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
            m2 <- apply(Y * Y * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
	} else if(family == "scaled-t"){
            ## scaled-t family
            m1 <- apply(abs(Y) * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
            m2 <- apply(Y * Y * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
        } else if(family == "gamma-binom") {
            ## gamma-binomial distribution
            ratio <- Y / X
            ratio[X == 0] <- mean(Y[X > 0] / X[X > 0])
            m1 <- apply(ratio * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
            m2 <- apply(ratio * ratio * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
        } else { 
  	    ## binomial
	    m1 <- apply(Y * ProbMat.full[idx, ], 1, sum) / apply(X * ProbMat.full[idx, ], 1, sum)
	    m2 <- NULL
	}
	MomentEstimate()
    }
    
    ## order the means
    od <-  apply(Mu, 1, order)
    Mu <- matrix(Mu[cbind(rep(seq_len(N), each = M), c(od))], ncol = M, byrow = TRUE)
    Sigma <- matrix(Sigma[cbind(rep(seq_len(N), each = M), c(od))], ncol = M, byrow = TRUE)

    assign("Mu", Mu, envir = parent.frame())
    assign("Sigma", Sigma, envir = parent.frame())
}

UpdateStates <- function() {
    Inherit()
    F1  <- matrix(0, nrow = K * S, ncol = I)
    F1.full <- matrix(0, nrow = N * M, ncol = I)
    Pi.full <- crossprod(t(designMap), Pi) ## N by S
    Pi.full <- tcrossprod(Pi.full, stateMap) ## N by M
    Pi.full <- Pi.full * V

    ## joint likelihood for each replicate and component
    for(m in seq(M)) {
      idx <- (m - 1) * N + seq_len(N)
      F1.full[idx, ] <- logdensity(Y, Mu[, m], Sigma[, m], X, family) + log(Pi.full[, m])
    }
    F1.full <- trimLogValue(F1.full)
    ##Update ProbMat.full
    totalF.full <- matrix(0, nrow = N, ncol = I)
    for(m in seq(M)) {
      idx <- seq(N) + (m - 1) * N
      F1.full[idx, ] <- exp(F1.full[idx, ])
      totalF.full <- totalF.full + F1.full[idx, ]
    }
    parlik <- sum(log(totalF.full))
    write.out(out, paste("Likelihood for component estimation: ", round(parlik, 3), sep = ""))
    totalF.full <- t(matrix(rep(c(t(totalF.full)), M), nrow = I))
    ProbMat.full <- F1.full / totalF.full
    ProbMat.full <- trimProbValue(ProbMat.full)

 
    ## compute F1
    ## recompute replicate density
    for(m in seq(M)) {
      idx <- (m - 1) * N + seq_len(N)
      F1.full[idx, ] <- exp(logdensity(Y, Mu[, m], Sigma[, m], X, family)) * V[, m]
    }
    ## convert to (NS) x I
    F1.full <- log(crossprod(unitMap, F1.full))
    for(s in seq(S)) {
      F1[(s - 1) * K + seq_len(K), ] <- crossprod(designMap, F1.full[(s - 1) * N + seq(N), ]) + log(Pi[, s])
    }
    F1 <- trimLogValue(F1)
    totalF <- matrix(0, nrow = K, ncol = I)
    for(s in seq_len(S)) {
        idx <- seq_len(K) + (s-1) * K
        F1[idx,] <- exp(F1[idx,])
        totalF <- totalF + F1[idx,]
    }
    totalF <- t(matrix(rep(c(t(totalF)), S), nrow = I))
    ProbMat <- F1 / totalF
    ProbMat <- trimProbValue(ProbMat)

    ## update Pi
    for(s in seq_len(S)) {
      idx <- seq_len(K) + (s-1) * K
      Pi[,s] <- apply(ProbMat[idx,], 1 ,mean)
    }

    ## update V
    tmp <- matrix(0, nrow = M * K, ncol = I)
    for(m in seq(M)) {
      tmp[(m - 1) * K + seq(K), ] <- ProbMat[(statemap[m] - 1) * K + seq(K), ]
    }
    EV <- matrix(0, nrow = N * M, ncol = I)
    for(m in seq(M)) {
      id1 <- (m - 1) * N + seq(N)
      id2 <- (m - 1) * K + seq(K)
      EV[id1, ] <- designMap %*% tmp[id2, ]
    }
    EV <- EV * c(V)
    EV <- EV + ProbMat.full
    
    V.new <- matrix(apply(EV, 1, sum), nrow = N)
    V.new <- V.new / tcrossprod(V.new %*% stateMap, stateMap)

    assign("V", V.new, envir = parent.frame())
    assign("Pi", Pi, envir = parent.frame())
    assign("ProbMat", ProbMat, envir = parent.frame())
    assign("ProbMat.full", ProbMat.full, envir = parent.frame())
}

PrintUpdate <- function() {
  Inherit()
  allerr <- sqrt(sum((ProbMat - ProbMat.true) ^ 2) / I / K / (S - 1))
  ## compute misclassification rate
  W.f <- matrix(0, nrow = K * S, ncol = J)
  for(s in seq_len(S))
    W.f[s + S * seq(0, K - 1),] <- W[seq_len(K) + K * (s - 1),]
  
  mc <- matchCluster(W.f, para$W, Zcond, para$Z, b.prob, para$non.id)
  
  write.out(out, paste("mis-class rate ", mc$mcr))
  write.out(out, paste("Error for W ",  round(mc$W.err, 3)))
  allmisclass <- c(allmisclass, mc$mcr)
  W.err <- c(W.err, mc$W.err)
  allari <- c(allari, mc$ari)
  Mu.err <- numeric(0)
  if(prod(dim(Mu) == dim(para$Mu)) == 1) {
    Mu.err <- sqrt(mean((Mu - para$Mu) ^ 2))
  }
  write.out(out, paste("ARI ", mc$ari))
  write.out(out, paste("loglik", totallik, "err", round(allerr, 3)))
  write.out(out, paste("Error for Mu", round(Mu.err, 3)))
  for(v in c("allerr", "allari", "allmisclass", "W.err", "Mu.err")) {
    assign(v, get(v), envir = parent.frame())
  }
}

Inherit <- function() {
    for(v in ls(envir = parent.frame(2))) {
        ## do not use get() since it will cause error for non-defined variables
        assign(v, parent.frame(2)[[v]], envir = parent.frame())
    }
}

logdensity <- function(y, mu, sigma, x = NULL, family) {
  if(family == "lognormal") {
    y <- log(y + 1)
    return(-(y - mu) ^ 2 / sigma / 2 - log(sigma) / 2 - log(2 * pi) / 2 )
  } else if(family == "scaled-t") {
    return(dt(y / mu, df = sigma, log = TRUE))
  } else if(family == "negbin") {
    return(dnbinom(y, mu = mu, size = sigma, log = TRUE))
  } else if(family == "gamma-binom") {
      a <- mu / (1 - mu) * sigma
      b <- sigma
      return(
          log(beta(a + y, x - y + b)) -
              log(beta(a, b)) + log(choose(x, y))
      )
  } else {
      return(dbinom(y, prob = mu, size = x, log = TRUE))
  }
}

MomentEstimate <- function() {
    for(v in c("m1", "m2", "Mu", "Sigma", "family", "m")) {
        assign(v, parent.frame()[[v]])
    }
    if(family == "lognormal") {
        Mu[, m] <- m1
        m2 <- m2 - m1 * m1
        m2[m2 < 0.01] <- 0.01
        Sigma[, m] <- m2
    } else if(family == "scaled-t") {
        func <- function(df) {
	    g1 <- gamma((df - 1) / 2)
	    g2 <- gamma(df / 2)
	    (df - 2) * g1 * g1 / g2 / g2 / pi
	}
	solve <- function(v) {
	    lower <- 2.01
	    upper <- 100
	    fun.up <- func(upper)
	    fun.low <- func(lower)
	    if(v >= fun.up) {
	        return(upper)
	    } else if(v <= fun.low) {
	        return(lower)
	    } else {
	        while(abs(fun.up - fun.low) > 0.01 & abs(upper - lower) > 0.01) {
		    med <- (upper + lower) / 2
		    fun.med <- func(med)
		    if(fun.med > v + 0.005) {
		        upper <- med
			fun.up <- fun.med
		    } else if(fun.med < v - 0.005){
		        lower <- med
			fun.low <- fun.med
		    } else {
		        return(med)
		    }
		}
		return((upper + lower) / 2)
	    }
	}
	## df
	Sigma[, m] <- sapply(m1 * m1 / m2, solve)
	## scale
	Mu[, m] <- sqrt(m2 * (1 - 2 / Sigma[, m]))
     } else if(family == "negbin") {
        Mu[, m] <- m1
        m2 <- m2 - m1 * m1
        m2 <- m1 / (m2 / m1 - 1)
        m2[m2 < 0] <- 100
        Sigma[, m] <- m2
    } else if(family == "gamma-binom") {
        m2 <- m2 - m1 * m1
        a.plus.b <- 1 / (m1 * (1 - m1)) - 1
        a <- m1 * a.plus.b
        b <- a.plus.b - a
        a[a < 0.01] <- 0.01
        b[b < 0.01] <- 0.01
        Mu[, m] <- a / (a + b)
        Sigma[, m] <- b
    } else {
        Mu[, m] <- m1
	Sigma[, m] <- 0
    }
    assign("Mu", Mu, envir = parent.frame())
    assign("Sigma", Sigma, envir = parent.frame())
}
