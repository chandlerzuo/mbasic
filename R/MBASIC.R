#' @name MBASIC
#' @title Bayesian clustering model for a state-space matrix.
#'
#' @description This function is designed to analyze general state-space models. The data consists of observations over I units under N experiments with K different conditions. There are S states for each experiment and unit.
#' @param Y An N by I matrix containing the data from N experiments across I observation units.
#' @param Gamma The data for background information. Default: \code{NULL}. See details for more information.
#' @param S An integer for the number of states.
#' @param statemap A vector the same length as the number of mixture components, and taking values from 1 to S representing the states of each component. Default: NULL. See details for more information.
#' @param fac A vector of levels repr1esenting the conditions of each replicate.
#' @param struct A K by J matrix indicating the structures of each cluster.
#' @param J The number of clusters to be identified.
#' @param family The distribution of family to be used. Either "lognormal", "negbin", "binom", "gamma-binom" or "scaled-t". See details for more information.
#' @param method A string for the fitting method, 'MBASIC' (default), 'PE-MC', 'SE-HC',  or 'SE-MC'. See details for more information.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for the relative increment in log-likelihood value in checking the algorithm convergence. Default: 1e-10.
#' @param tol.par Tolerance for the relative error in parameter updates in checking the algorithm convergence. Default: 1e-5.
#' @param zeta The initial value for the proportion of units that are not clustered. Default: 0.1. If 0, no singleton cluster is fitted.
#' @param out The file directory for writing fitting information in each E-M iteration. Default: NULL (no information is outputted).
#' @param verbose A boolean variable indicating whether intermediate model fitting metrics should be printed. Default: FALSE.
#' @param min.count The minimum count threshold for each component. This argument is only used when \code{family='negbin'}. If it is a single number, it is the threshold for all components >= 2. If it is a vector, it is the threshold for each component.
#' @param para A list object that contains the true model parameters. Default: NULL. See details for more information.
#' @param initial Either a \code{\link{list}} or \linkS4class{MBASICFit} object that provides initial values for model parameters. Default: NULL.
#' @details
#' MBASIC assumes that there are S underlying states for each expeirment and each loci. A single state may also include multiple mixture components, indexed by m. In total, we can have M mixture components. The mapping from mixture components to the states are provided by \code{statemap}. By default, \code{statemap=NULL}, in which case each state has only one component, and M=S.\cr
#' Function MBASIC currently supports five different distributional families: log-normal, negative binomial, binomial, gamma-binomial and scaled-t distributions. This should be specified by the \code{family} argument.\cr
#' For the log-normal distributions, log(Y+1) is modeled as normal distributions. For experiment n, if locus i has component m, distribution for log(Y[n,i]+1) is N(Mu[n,m]*Gamma[n,i+I(m-1)], Sigma[n,m]).\cr
#' For the negative binomial distributions, the meanings of Mu and Sigma are different. For experiment n, if locus i has component m, distribution of Y[n,i]-min.count[m] is NB(Mu[n,m]*Gamma[n,i+I(m-1)], Sigma[n,m]). In this package, NB(mu, a) denotes the negative-binomial distribution with mean mu and size a (i.e. the variance is mu*(1+mu/a)). Notice that if a single value of 'min.count' is provided, it will be converted to a vector of \code{c(0, rep(min.count, M-1))}.\cr
#' For the binomial distribution, for experiment n, if locus i has component m, distribution for Y[n,i] is Binom(Gamma[n,i], Mu[n,m]).\cr
#' For the gamma-binomial distribution, for experiment n, if locus i has component m, distribution for Y[n,i] is Binom(Gamma[n,i], p) where p follows a gamma prior of gamma(Mu[n,m], Sigma[n,m]).\cr
#' For the scaled-t distribution, for experiment n, if locus i has component m, distribution for Y[n,i]/Gamma[n,i+I(m-1)]/Mu[n,m] is t distribution with Sigma[n,m] degrees of freedom.\cr
#' The \code{Gamma} parameter encodes the background information for all N experiments, I units and M components. It can be a matrix with dimension K by I * M, where the background datum for experiment n, unit i and component m is Gamma[n,i+I*(m-1)]. If in the input \code{Gamma=NULL}, then it is regenerated as a matrix of entries 1 with dimension N x IM. If in the input \code{Gamma} is a N x I matrix, then this function adds I(M-1) columns of all 1s to this matrix.\cr
#'  The \code{method} argument determines what fitting method will be used. The default is 'MBASIC', where the states and the clustering are simultaneously estimated. 'SE-HC' and 'SE-MC' methods use 2-step algorithms. In the first step, both estimate the states for each unit by an E-M algorithm for each experiment. In the second step, 'SE-HC' uses hierarchical clustering to cluster the units, while 'SE-MC' uses function \code{\link{MBASIC.state}} to identify clusters.\cr
#' The \code{para} argument takes a list object that is supposed to include the following fields:
#'\tabular{ll}{
#' W \tab A K by (J*S) matrix. The (k,J*(s-1)+j)-th entry is the probability that the units in cluster j has state s in the k-th experiment.\cr
#' Z \tab An I by J matrix. The (i,j)-th entry is the indicator whether the i-th unit belongs to cluster j.\cr
#' Theta \tab A K by (I*S) matrix. The (k,I*(s-1)+i)-th entry is the probability that the i-th unit has state s in the k-th experiment.\cr
#' non.id \tab A binary vector of length I. The i-th entry is the indicator whether the i-th unit does not belong to any cluster.
#' }
#' This argument is intended to carry the true parameters in simulation studies. If it is not null, then the model also computes a number of metrics that describes the error in model fitting. Users should be cautious that the order of the rows and columns of matrices in the fields of para should match the Y matrix.
#' @return An object of class \linkS4class{MBASICFit}.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' ## Simulate a dataset
#' dat.sim <- MBASIC.sim(xi = 2, family = "lognormal", I = 1000, fac = rep(1:10, each = 2), J = 3, S = 3, zeta = 0.1)
#' ## Fit the model
#' dat.sim.fit <- MBASIC(Y = dat.sim$Y, S = 3, fac = rep(1:10, each = 2), J = 3, maxitr = 3, para = NULL, family = "lognormal", method = "MBASIC", zeta = 0.1, tol = 1e-6)
#' @useDynLib MBASIC
#' @export
MBASIC <- function(Y, Gamma = NULL, S, fac, J=NULL, maxitr = 100, struct = NULL, para = NULL,  family="lognormal", method = "MBASIC", zeta = 0.1, min.count = 0, tol = 1e-10, tol.par = 0.001, out = NULL, verbose = FALSE, statemap = NULL, initial = NULL) {

  Mu.init <- Sigma.init <- V.init <- ProbMat.init <- W.init <- Z.init <- b.init <- P.init <- NULL

  if(!is.null(initial)) {
    if(class(initial) == "list") {
      for(var.init in c("Mu", "Sigma", "V", "ProbMat", "W", "Z", "b", "P")) {
        if(var.init %in% names(initial)) {
          assign(paste(var.init, "init", sep = "."), initial[[var.init]])
        }
      }
    } else if(class(initial) == "MBASICFit") {
      for(var.init in c("Mu", "Sigma", "V", "W", "Z", "b", "P")) {
        if(var.init %in% names(initial)) {
          if(length(slot(initial, var.init)) > 0) {
            assign(paste(var.init, "init", sep = "."), slot(initial, var.init))
          }
        }
      }
      ProbMat.init <- initial@Theta
    } else {
      warning("'initial' is not list or class 'MBASICFit'. Its value is skipped.")
    }
  }
  
  if(verbose) {
    write.out(out, "Started")
  }
  if(! method %in% c("SE-HC", "SE-MC", "PE-MC", "MBASIC")) {
    stop("Error: 'method' must be one of 'SE-HC', 'SE-MC', 'PE-MC' or 'MBASIC'.")
  }

  if(! family %in% c("lognormal", "negbin", "binom", "gamma-binom", "scaled-t")) {
    stop("Error: 'family' must be one of 'lognormal', 'negbin', 'binom', 'gamma-binom', 'scaled-t'.")
  }

  if(is.null(statemap)) {
    ## states and components have one-to-one matching
    statemap <- seq(S)
  }

  M <- length(statemap)
  stateMap <- matrix(0, nrow = M, ncol = S)
  stateMap[cbind(seq(M), statemap)] <- 1

  if(length(min.count) == 1) {
    min.count <- rep(min.count, M)
    min.count[statemap == 1] <- 0
  } else if(is.null(min.count)) {
    min.count <- rep(0, M)
  } else if(length(min.count) != M) {
    warning("Length of 'min.count' is not the same as 'S' or 'stateMap'. No minimum count threshold value is set.")
    min.count <- rep(0, M)
  }
  if(family != "negbin" & sum(min.count != 0) > 0) {
    warning("'min.count' is only supported for the negative binomial distribution. No minimum count threshold value is set.")
    min.count <- rep(0, M)
  }
  if(sum(min.count < 0) > 0) {
    warning("Negative values in 'min.count' are reset as 0.")
    min.count[min.count < 0] <- 0
  }
  
  ## prespecified
  K <- length(unique(fac))
  facNames <- as.character(unique(fac))
  facMap <- seq(K)
  names(facMap) <- facNames
  fac <- as.character(fac)
  fac <- facMap[as.character(fac)]
  I <- ncol(Y)
  N <- nrow(Y)
  if(length(fac) != N)
    stop("Error: total number of replicates do not match with the number of rows in Y")

  if(family == "binom") {
    if(is.null(Gamma)) {
      Gamma <- matrix(apply(Y, 1, max), nrow = N, ncol = I)
    } else if(nrow(Gamma) != N | ncol(Gamma) != I) {
      stop("Dimension of 'Gamma' must be the same as 'Y'.")
    } else {
      if(sum(Gamma < Y) > 0) {
        warning("Some 'Gamma' is smaller than 'Y'. Those covariates are automatically augmented.")
        Gamma[Gamma < Y] <- Y[Gamma < Y]
      }
    }
    Gamma[Gamma == 0] <- 1
  }
  
  if(is.null(Gamma)) {
    Gamma <- matrix(0.5, nrow = N, ncol = I * M)
  } else if(nrow(Gamma) == N & ncol(Gamma) == I) {
    Gamma.add <- matrix(apply(Gamma, 1, mean), nrow = N, ncol = I * (M - 1))
    Gamma <- cbind(Gamma, Gamma.add)
  } else if(nrow(Gamma) != N | ncol(Gamma) != I * M) {
    stop("Error: structure of 'Gamma' is not correct. See details.")
  }
  if(min(Gamma) == 0) {
    Gamma <- Gamma + min(Gamma[Gamma > 0])
  }
  
  ## scale the Gamma matrix
  scaleMat <- matrix(1, nrow = N, ncol = M)
  for(m in seq(M)) {
    scaleMat[, m] <- apply(Gamma[, I * (m - 1) + seq(I)], 1, mean)
    Gamma[, I * (m - 1) + seq(I)] <- Gamma[, I * (m - 1) + seq(I)] / scaleMat[, m]
  }
  
  if(is.null(struct)) {
    if(is.null(J))
      stop("Error: either struct or J must not be missing.")
    struct <- matrix(seq_len(K), nrow = K, ncol = J)
  } else {
    if(is.null(J))
      J <- ncol(struct)
    J <- sum(J)
    if(ncol(struct)!= sum(J) | nrow(struct) != K)
      stop("Error: the dimension of struct is inconsistent with grouping structure!")
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
    J +                   ## probz and zeta
      sum(apply(struct, 2, function(x) length(unique(x)))) * (S - 1) + ## W
        N * M +       ## Mu
          N * (M - S) ## V

  if(family != "binom") {
    numpar <- numpar + N * M ## Sigma
  }
  
  outitr <- 0
  totallik <- oldlik <- -Inf
  alllik <- allerr <- allzeta <- allmisclass <- matchId1 <- W.err <- matchId2 <- allari <- numeric(0)
  maxlik <- -Inf

  if(verbose) {
    write.out(out, "Initialize parameters")
  }
  
  ## initialize distributions
  V <- Sigma <- Mu <- matrix(0, nrow = N, ncol = M)

  if(is.null(Mu.init) | is.null(Sigma.init)) {
    InitDist()
  } else {
    Mu <- Mu.init
    Sigma <- Sigma.init
  }
  
  b <- rep(0, I)
  
  ## initialize the matrices by K-means
  ## in constructing Z, cluster all locis
  ProbMat <- matrix(0, nrow = K * S, ncol = I)
  if(M == S) {
    V.init <- matrix(1, ncol = M, nrow = N)
  }
  if(is.null(ProbMat.init) | is.null(V.init)) {
    InitStates()
  } else {
    ProbMat <- ProbMat.init
    V <- V.init
    Pi <- matrix(0, nrow = K, ncol = S)
    for(s in seq_len(S)) {
      idx <- seq_len(K) + (s-1) * K
      if(length(idx) > 1) {
        Pi[,s] <- apply(ProbMat[idx,], 1 ,mean)
      } else {
        Pi[, s] <- mean(ProbMat[idx, ])
      }
    }
    ProbMat.full <- matrix(0, nrow = M * N, ncol = I)
    for(m in seq(M)) {
      ProbMat.full[(m - 1) * N + seq(N), ] <- V[, m] * (designMap %*% ProbMat[(statemap[m] - 1) * K + seq(K), ])
    }
  }

  if(!is.null(para)) {
    ProbMat.true <- matrix(0, nrow = S * K, ncol = I)
    for(s in seq(S)) {
      ProbMat.true[(s - 1) * K + seq(K), ] <- as.numeric(para$Theta == s)
    }
    W.true <- para$W
    Z.true <- para$Z
    nonid.true <- para$non.id
    Mu.true <- para$Mu
  }

  if(method != "MBASIC") {
    ## SE-HC, PE-MC or SE-MC method
    
    ## em step to estimate Theta
    allpar <- c(c(V), c(Mu), c(Sigma), c(Pi))
    oldpar <- 0
    alllik <- numeric(0)
    totallik <- oldlik <- maxlik <- -Inf
    conv <- FALSE
    for(itr in seq(maxitr)) {
      oldlik <- totallik
      
      ## M step
      UpdateDist()
      
      ## E step
      UpdateStates()
      oldpar <- allpar
      allpar <- c(c(V), c(Mu), c(Sigma), c(Pi))
      alllik <- c(alllik, totallik)
      if(maxlik < totallik) {
        assignBest(c("ProbMat", "Mu", "Sigma", "Pi", "V"))
        maxlik <- totallik
      }

      ## check for convergence
      if(max(na.omit(abs(oldpar / allpar - 1))) < tol.par) {
        conv <- TRUE
        break
      }
      if(itr > 10 & oldlik < totallik & totallik - oldlik < tol * abs(oldlik)) {
        conv <- TRUE
        break
      }
      if(itr > 100 & which.max(alllik) < itr / 2) {
        conv <- TRUE
        break
      }
    } ## finish iteration
    getBest(c("ProbMat", "Mu", "Sigma", "Pi", "V"))
    
    Theta <- matrix(-1, nrow = K, ncol = I)
    for(k in seq_len(K)) {
      idx <- k + K * (seq_len(S) - 1)
      Theta[k,] <- apply(ProbMat[idx,], 2, which.max)
    }

    if(!is.null(para))
      allerr <- sqrt(sum((ProbMat - ProbMat.true) ^ 2) / I / K / S)
    
    if(method != "PE-MC") {
      ret <- MBASIC.state(Theta, J=J, zeta = zeta, struct = struct, method = method, maxitr = maxitr, tol = tol, tol.par = tol.par, para = para, out = out, W.init = W.init, Z.init = Z.init, P.init = P.init, b.init = b.init)
      
      conv <- (conv & ret@converged)
      
      ## Pi is the proportion for components in the k experiment to have state s
      ## Pi is different from Z. Z is the posterior probability.

      Mu.err <- numeric(0)
      if(prod(dim(Mu) == dim(Mu.true)) == 1) {
        Mu.err <- sqrt(mean((Mu - Mu.true) ^ 2))
      }

      write.out(out, paste("mis-class rate ", ret@MisClassRate))
      write.out(out, paste("Error for W ",  round(ret@W.err, 3)))
      write.out(out, paste("ARI ", ret@ARI))
      write.out(out, paste("loglik", round(tail(ret@alllik, 1), 3), "err", round(allerr, 3)))
      write.out(out, paste("Error for Mu", round(Mu.err, 3)))

      W <- ret@W
      rownames(ProbMat) <- facNames
      rownames(W) <- rep(facNames, S)
      rownames(Mu) <- rownames(Sigma) <- rownames(V) <- facNames[fac]

      return(new("MBASICFit",
                 Theta = ProbMat,
                 W = W,
                 Z = ret@Z,
                 V = V,
                 b = ret@b,
                 clustProb = ret@clustProb,
                 lik = ret@lik,
                 alllik = ret@alllik,
                 zeta = ret@zeta,
                 Mu = Mu / scaleMat,
                 Sigma = Sigma,
                 probz = ret@probz,
                 P = ret@P,
                 converged = conv,
                 Theta.err = allerr,
                 ARI = ret@ARI,
                 W.err = ret@W.err,
                 MisClassRate = ret@MisClassRate,
                 Mu.err = Mu.err,
                 Iter = ret@Iter
                 )
             )
    }
  }

  ## initialize W, Z, b
  ## ProbMat <- D.rep %*% ProbMat
  if(is.null(W.init) | is.null(Z.init)) {
    InitWZb()
  } else {
    W <- W.init
    Z <- Z.init
  }
  if(!is.null(b.init)) {
  }
  predZ <- Zcond <- Z
  
  ## initialize p, probz
  if(is.null(P.init)) {
    P <- matrix(0, nrow = I, ncol = S)
    for(s in seq_len(S)) {
      idx <- seq_len(K) + K * (s - 1)
      P[, s] <- apply(ProbMat[idx, ], 2, mean)
    }
  } else {
    P <- P.init
  }
  probz <- apply(rbind(Z, diag(rep(1, J))), 2, mean)
  
  ## EM algorithm
  ## change storage modes for C modules
  storage.mode(statemap) <- "integer"
  storage.mode(fac) <- "integer"
  oldpar <- 0
  newpar <- c(c(W), probz, zeta, c(P), c(V), c(Mu), c(Sigma))
  alllik <- numeric(0)
  totallik <- oldlik <- maxlik <- -Inf
  conv <- FALSE

  for(outitr in seq_len(maxitr)) {
    if(outitr == 1 | method == "MBASIC") {
      ## only compute the PDF once if the method is PE-MC
      PDF <- matrix(0, nrow = N * M, ncol = I)
      for(m in seq_len(M)) {
        Gamma.m <- Gamma[, (m - 1) * I + seq(I)]
        PDF[seq(N) + (m - 1) * N, ] <- logdensity(Y, Mu[, m], Sigma[, m], Gamma.m, family, min.count[m])
      }
      PDF <- trimLogValue(PDF)
    }
   
    Theta <- matrix(-1, nrow = K, ncol = I)
    for(k in seq_len(K)) {
      idx <- k + K * (seq_len(S) - 1)
      Theta[k,] <- apply(ProbMat[idx,], 2, which.max)
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
    P <- estep.result[["P"]]

    clustOrder <- .orderCluster(W, struct)
    W <- W[, clustOrder, drop = FALSE]
    W <- .structure(W, struct)
    probz <- probz[clustOrder]
    Zcond <- Zcond[, clustOrder, drop = FALSE]
    predZ <- predZ[, clustOrder, drop = FALSE]
    
    if(method != "PE-MC") {
      ## skip this step if the method is PE-MC
      ## M-step for Mu and Sigma
      UpdateDist()
    }
    
    oldlik <- totallik
    totallik <- .Call("loglik", W, P, V, zeta, probz, PDF, fac - 1, statemap - 1, package="MBASIC")
    if(verbose) {
      write.out(out, paste("itr", outitr, "lik", round(tail(totallik, 1), 2), "zeta", round(zeta, 2)))
    }
    alllik <- c(alllik, totallik)
    allzeta <- c(allzeta, zeta)
    if(length(para) > 1 & verbose) {
      PrintUpdate()
    }
    
    if(maxlik < totallik) {
      maxlik <- totallik
      assignBest(c("ProbMat", "Theta", "W", "V", "P", "Mu", "Sigma", "zeta", "probz", "predZ", "Zcond", "b.prob"))
    }
    oldpar <- newpar
    newpar <- c(c(W), probz, zeta, c(P), c(V), c(Mu), c(Sigma))
    if(max(na.omit(abs(newpar / oldpar - 1))) < tol.par) {
      conv <- TRUE
      break
    }
    if(outitr > 10 & oldlik < totallik & totallik - oldlik < tol * abs(oldlik)) {
      conv <- TRUE
      break
    }
    if(outitr > 100 & which.max(alllik) < outitr  - 100) {
      conv <- TRUE
      break
    }
    
  } ## finish outer loop
  
  getBest(c("ProbMat", "Theta", "W", "V", "P", "Mu", "Sigma", "zeta", "probz", "predZ", "Zcond", "b.prob"))

  Mu.err <- W.err <- numeric(0)
  if(length(para) > 1)
    PrintUpdate()

  rownames(ProbMat) <- rep(facNames, S)
  rownames(Mu) <- rownames(Sigma) <- rownames(V) <- facNames[fac]
  rownames(W) <- rep(facNames, S)
  
  new("MBASICFit",
      Theta = ProbMat,
      W = W,
      V = V,
      Z = Z,
      b = b.prob,
      clustProb = cbind(b.prob, Zcond * (1 - b.prob)),
      aic = - 2 * tail(maxlik, 1) + 2 * numpar,
      bic = - 2 * tail(maxlik, 1) + log(N * I) * numpar,
      aicc = -2 * tail(maxlik, 1) + 2 * numpar + 2 * numpar * (numpar + 1) / (N * I - numpar - 1),
      alllik = alllik,
      lik = maxlik,
      zeta = zeta,
      Mu = Mu / scaleMat,
      Sigma = Sigma,
      probz = probz,
      P = P,
      converged = conv,
      Theta.err = allerr,
      ARI = tail(allari, 1),
      W.err = tail(W.err, 1),
      MisClassRate = tail(allmisclass, 1),
      Struct = struct,
      Mu.err = Mu.err,
      Iter = outitr
      )
}

InitStates <- function() {
  Inherit(c("S", "statemap", "V", "K", "M", "I", "N", "Y", "Sigma", "Mu", "family", "unitMap", "designMap", "stateMap", "Gamma", "min.count"))
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
    Gamma.m <- Gamma[, I * (m - 1) + seq(I)]
    idx <- (m - 1) * N + seq_len(N)
    F1.full[idx,] <- logdensity(Y, Mu[, m], Sigma[, m], Gamma.m, family, min.count[m])
  }
  F1.full <- trimLogValue(F1.full)
  F1.full <- exp(F1.full)
  ## convert into (NS) x I
  F1.tmp <- log(crossprod(unitMap, F1.full))
  
  ## initialize ProbMat.full
  for(m in seq(M)) {
    idx <- (m - 1) * N + seq(N)
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

  Return(c("ProbMat.full", "ProbMat", "V", "Pi"))
}

InitDist <- function() {
  Inherit(c("M", "Y", "family", "Mu", "Sigma", "Gamma", "min.count"))
  Gamma.m <- Gamma[, seq(ncol(Y))]
  
  for(m in seq(M)) {
    m3 <- NULL
    thr <- min.count[m]
    if(family == "lognormal") {
      Y.sec <- c(Y)[c(Y) <= quantile(c(Y), m / M) & c(Y) >= quantile(c(Y), (m - 1) / M)]
      Y.sec <- log(Y.sec + 1)
    } else if(family  == "negbin") {
      Y.sec <- c(Y)[c(Y) < quantile(c(Y), m / M) & c(Y) > quantile(c(Y), (m - 1) / M)]
      Y.sec <- Y.sec - thr
      Gamma.sec <- c(Gamma * Gamma)[c(Y) < quantile(c(Y), m / M) & c(Y) > quantile(c(Y), (m - 1) / M)]
      m3 <- mean(Gamma.sec)
    } else if(family == "scaled-t") {
      Y.sec <- abs(Y)
      Y.sec <- c(Y.sec)[c(Y.sec) < quantile(c(Y.sec), m / M) & c(Y.sec) > quantile(c(Y.sec), (m - 1) / M)]
    } else if(family == "gamma-binom") {
      ## gamma-binomial distribution
      ratio <- Y / Gamma.m
      ratio[Gamma.m == 0] <- mean(Y[Gamma.m > 0] / Gamma.m[Gamma.m > 0])
      Y.sec <- c(ratio)[ratio <= quantile(ratio, m / M) & ratio >= quantile(ratio, (m - 1) / M)]
    } else {
      ratio <- Y / Gamma.m
      id1 <- which(Gamma.m > 0)
      ratio <- ratio[id1]
      id <- which(ratio <= quantile(ratio, m / M) & ratio >= quantile(ratio, (m - 1) / M))
      Y.sec <- sum(Y[id1[id]]) / sum(Gamma.m[id1[id]])
    }
    m1 <- mean(Y.sec)
    m2 <- var(Y.sec)
      
    MomentEstimate()
  }
  Return(c("Mu", "Sigma"))
}

UpdateDist <- function() {
  Inherit(c("M", "N", "family", "Y", "ProbMat.full", "Mu", "Sigma", "Gamma", "min.count"))

  I <- ncol(Y)
  for(m in seq_len(M)) {
    idx <- seq(N) + (m - 1) * N
    Gamma.m <- Gamma[, I * (m - 1) + seq(I)]
    thr <- min.count[m]
    m3 <- NULL
    if(family == "lognormal") {
      m1 <- apply(log(Y + 1) * ProbMat.full[idx, ], 1, sum) / apply(Gamma.m * ProbMat.full[idx, ], 1, sum)
      m2 <- apply((log(Y + 1) - m1 * Gamma.m) * (log(Y + 1) - m1 * Gamma.m) * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
    } else if(family == "negbin"){
      ## negative binomial family
      m1 <- apply((Y - thr) * ProbMat.full[idx, ], 1, sum) / apply(Gamma.m * ProbMat.full[idx, ], 1, sum)
      m2 <- apply((Y - thr) * (Y - thr) * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
      m3 <- apply(Gamma.m * Gamma.m * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
    } else if(family == "scaled-t"){
      ## scaled-t family
      m1 <- apply(abs(Y / Gamma.m) * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
      m2 <- apply(Y * Y / Gamma.m / Gamma.m * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
    } else if(family == "gamma-binom") {
      ## gamma-binomial distribution
      ratio <- Y / Gamma.m
      ratio[Gamma.m == 0] <- mean(Y[Gamma.m > 0] / Gamma.m[Gamma.m > 0])
      m1 <- apply(ratio * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
      m2 <- apply(ratio * ratio * ProbMat.full[idx, ], 1, sum) / apply(ProbMat.full[idx, ], 1, sum)
    } else { 
      ## binomial
      m1 <- apply(Y * ProbMat.full[idx, ], 1, sum) / apply(Gamma.m * ProbMat.full[idx, ], 1, sum)
      m2 <- NULL
    }
    MomentEstimate()
  }
  
  ## order the means
  Mu <- Mu + rep(min.count, each = nrow(Mu))
  od <-  apply(Mu, 1, order)
  Mu <- matrix(Mu[cbind(rep(seq_len(N), each = M), c(od))], ncol = M, byrow = TRUE)
  Mu <- Mu - rep(min.count, each = nrow(Mu))
  Sigma <- matrix(Sigma[cbind(rep(seq_len(N), each = M), c(od))], ncol = M, byrow = TRUE)

  Return(c("Mu", "Sigma"))
}

UpdateStates <- function() {
  Inherit(c("K", "S", "I", "M", "N", "Y", "Gamma", "Mu", "Sigma", "family", "V", "unitMap", "designMap", "stateMap", "Pi", "statemap", "out", "verbose", "min.count"))
  F1  <- matrix(0, nrow = K * S, ncol = I)
  F1.full <- matrix(0, nrow = N * M, ncol = I)

  ## compute F1
  ## recompute replicate density
  for(m in seq(M)) {
    idx <- (m - 1) * N + seq_len(N)
    Gamma.m <- Gamma[, seq(I) + I * (m - 1)]
    F1.full[idx, ] <- exp(logdensity(Y, Mu[, m], Sigma[, m], Gamma.m, family, min.count[m])) * V[, m]
  }
  F1.tmp <- F1.full
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

  ## joint likelihood for each replicate and component
  totallik <- sum(log(totalF))
  if(verbose)
    write.out(out, paste("Likelihood for component estimation: ", round(totallik, 3), sep = ""))

  ## ProbMat.full
  ProbMat.full <- F1.tmp
  for(s in seq(S)) {
    Vwei <- 0
    for(m in which(statemap == s)) {
      Vwei <- Vwei + ProbMat.full[(m - 1) * N + seq(N), ]
    }
    for(m in which(statemap == s)) {
      ProbMat.full[(m - 1) * N + seq(N), ] <- ProbMat.full[(m - 1) * N + seq(N), ] / Vwei
    }
  }
  ProbMat.aug <- matrix(0, nrow = N * M, ncol = I)
  for(m in seq(M)) {
    ProbMat.aug[(m - 1) * N + seq(N), ] <- designMap %*% ProbMat[(statemap[m] - 1) * K + seq(K), ]
  }
  ProbMat.aug <- trimProbValue(ProbMat.aug)
  ProbMat.full <- ProbMat.full * ProbMat.aug
  ProbMat.full <- trimProbValue(ProbMat.full)

  ## update Pi
  for(s in seq_len(S)) {
    idx <- seq_len(K) + (s-1) * K
    if(length(idx) > 1) {
      Pi[,s] <- apply(ProbMat[idx,], 1 ,mean)
    } else {
      Pi[, s] <- mean(ProbMat[idx, ])
    }
  }

  ## update V
  EV <- ProbMat.aug * c(V) + ProbMat.full
  V.new <- matrix(apply(EV, 1, sum), nrow = N)
  V.new <- V.new / tcrossprod(V.new %*% stateMap, stateMap)

  Return(c("V", "Pi", "ProbMat", "ProbMat.full", "totallik"))
}

PrintUpdate <- function() {
  Inherit(c("ProbMat", "ProbMat.true", "I", "K", "S", "J", "W", "W.true", "Z.true", "nonid.true", "b.prob", "Zcond", "allmisclass", "W.err", "allari", "Mu", "Mu.true", "totallik", "out", "verbose"))
  allerr <- sqrt(sum((ProbMat - ProbMat.true) ^ 2) / I / K / S)
  ## compute misclassification rate
  W.f <- matrix(0, nrow = K * S, ncol = J)
  for(s in seq_len(S))
    W.f[s + S * seq(0, K - 1),] <- W[seq_len(K) + K * (s - 1),]
  
  mc <- matchCluster(W.f, W.true, Zcond, Z.true, b.prob, nonid.true)
  
  allmisclass <- c(allmisclass, mc$mcr)
  W.err <- c(W.err, mc$W.err)
  allari <- c(allari, mc$ari)
  Mu.err <- numeric(0)
  if(prod(dim(Mu) == dim(Mu.true)) == 1) {
    Mu.err <- sqrt(mean((Mu - Mu.true) ^ 2))
  }
  if(verbose) {
    write.out(out, paste("mis-class rate ", mc$mcr))
    write.out(out, paste("Error for W ",  round(mc$W.err, 3)))
    write.out(out, paste("ARI ", mc$ari))
    write.out(out, paste("loglik", totallik, "err", round(allerr, 3)))
    write.out(out, paste("Error for Mu", round(Mu.err, 3)))
  }
  Return(c("allerr", "allari", "allmisclass", "W.err", "Mu.err"))
}

Inherit <- function(vnames = NULL) {
  if(is.null(vnames)) {
    vnames <- ls(envir = parent.frame(2))
  }
  for(v in vnames) {
    ## do not use get() since it will cause error for non-defined variables
    assign(v, parent.frame(2)[[v]], envir = parent.frame())
  }
}

Return <- function(vnames = NULL) {
  if(is.null(vnames)) {
    vnames <- ls(envir = parent.frame())
  }
  for(v in vnames) {
    ## do not use get() since it will cause error for non-defined variables
    assign(v, parent.frame()[[v]], envir = parent.frame(2))
  }
}

logdensity <- function(y, mu, sigma, gamma, family, thr) {
  if(family == "lognormal") {
    y <- log(y + 1)
    return(-(y - mu * gamma) ^ 2 / sigma / 2 - log(sigma) / 2 - log(2 * pi) / 2 )
  } else if(family == "scaled-t") {
    return(dt(y / mu / gamma, df = sigma, log = TRUE) - log(mu * gamma))
  } else if(family == "negbin") {
    return(dnbinom(y - thr, mu = mu * gamma, size = sigma, log = TRUE))
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
  Inherit(c("m1", "m2", "m3", "Mu", "Sigma", "family", "m", "thr"))
  if(family == "lognormal") {
    Mu[, m] <- m1
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
    Sigma[, m] <- sapply(m1 * m1 / (m2 + m1 * m1), solve)
    ## scale
    Mu[, m] <- sqrt((m2 + m1 * m1) * (Sigma[, m] - 2) / Sigma[, m])
  } else if(family == "negbin") {
    m1[m1 < 0.001] <- 0.001
    Mu[, m] <- m1
    m2 <- (m1 * m1 * m3) / (m2 - m1 * m1 - m1 * m1 * m3)
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
  Return(c("Mu", "Sigma"))
}

assignBest <- function(varnames) {
  for(v in varnames) {
    bestname <- paste("best", v, sep = "")
    assign(bestname, get(v, envir = parent.frame()), envir = parent.frame())
  }
}

getBest <- function(varnames) {
  for(v in varnames) {
    bestname <- paste("best", v, sep = "")
    assign(v, get(bestname, envir = parent.frame()), envir = parent.frame())
  }
}

InitWZb <- function() {
  Inherit(c("ProbMat", "zeta", "I", "J", "S", "K", "struct"))
  km.fit <- kmeans(t(ProbMat), centers = J)
  Z <- matrix(0, nrow = I, ncol = J)
  Z[cbind(seq(I), km.fit$cluster)] <- 1
  W <- t(km.fit$centers)
  d <- apply(ProbMat - tcrossprod(W, Z), 1, function(x) sum(x * x))
  thr <- quantile(d, 1 - zeta)
  id <- which(d < thr)
  b <- rep(1, I)
  b[id] <- 0
  clustOrder <- .orderCluster(W, struct)
  W <- W[, clustOrder, drop = FALSE]
  Z <- Z[, clustOrder, drop = FALSE]
  W <- .structure(W, struct)
  b.prob <- b

  Return(c("W", "Z", "b.prob"))

}

#' @name MBASIC.full
#' @title Simultaneously fitting MBASIC model for different numbers of clusters.
#' @description Simultaneously fitting MBASIC model for different numbers of clusters.
#' @param struct A list of K by J matrix indicating the structures of each cluster. Default: NULL.
#' @param J A vector of the numbers of clusters to be identified.
#' @param ncores The number of CPUs to be used in parallele. Default: 1.
#' @param out The prefix of the file names to write model output information. Default: NULL, when the information is written to the R output.
#' @param ... Other parameters that are passed to \code{\link{MBASIC}}.
#' @details
#' This function fits MBASIC models with different numbers of clusters and/or different structures in parallel.
#' @return A list object including the following fields:
#' \tabular{ll}{
#' allFits \tab A list of \linkS4class{MBASICFit} objects with different numbers of clusters and/or structures.\cr
#' BestFit \tab The best model with minimum BIC.\cr
#' Time \tab The time used for the model fitting.\cr}
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @seealso \code{\link{MBASIC}}
#' @examples
#' ## Simulate a dataset
#' dat.sim <- MBASIC.sim(xi = 2, family = "lognormal", I = 1000, fac = rep(1:10, each = 2), J = 3, S = 3, zeta = 0.1)
#' ## Fit the model
#' dat.sim.fit <- MBASIC.full(Y = dat.sim$Y, S = 3, fac = rep(1:10, each = 2), J = 3:5, maxitr = 3, para = NULL, family = "lognormal", method = "MBASIC", zeta = 0.1, tol = 1e-6, tol.par = 0.001)
#' @useDynLib MBASIC
#' @import doMC
#' @export
MBASIC.full <- function(J=NULL, struct = NULL, out = NULL, ncores = 1, ...) {
  t0 <- Sys.time()
  allJs <- J
  allstructs <- struct
  if(is.null(allstructs)) {
    nmodels <- length(allJs)
  } else {
    nmodels <- length(allstructs)
  }
  if(nmodels == 0) {
    stop("Error: at least J or struct must be not NULL.")
  }

  registerDoMC(min(c(ncores, nmodels)))
  allfits <- foreach(i = seq(nmodels)) %dopar% {
    if(!is.null(allstructs)) {
      struct <- allstructs[[i]]
    }
    if(!is.null(allJs)) {
      J <- allJs[i]
    }
    if(!is.null(out)) {
      out <- paste(out, "model", i, ".txt", sep = "")
    }
    MBASIC(J = J, struct = struct, out = out, ...)
  }
  bestBIC <- Inf
  bestFit <- NULL
  for(fit in allfits) {
    if(fit@bic < bestBIC) {
      bestBIC <- fit@bic
      bestFit <- fit
    }
  }
  return(list(allFits = allfits, BestFit = bestFit,
              Time = as.numeric(Sys.time() - t0, units = "secs")))
}

