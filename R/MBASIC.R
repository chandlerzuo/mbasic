#' @name MBASIC
#' @title Bayesian clustering model for a state-space matrix.
#'
#' @description This function is designed to analyze general state-space models. The data consists of observations over I units under N experiments with K different conditions. There are S states for each experiment and unit.
#' @param Y An N by I matrix containing the data from N experiments across I observation units.
#' @param S An integer for the number of states.
#' @param fac A vector of levels repr1esenting the conditions of each replicate.
#' @param struct A K by J matrix indicating the structures of each cluster.
#' @param J The number of clusters to be identified.
#' @param family The distribution of family to be used. Either "lognormal", "negbin" or "gamma-binom". See details for more information.
#' @param method A string for the fitting method, 'MBASIC' (default), 'PE-MC', 'SE-HC',  or 'SE-MC'. See details for more information.
#' @param para A list object that contains the true model parameters. Default: NULL. See details for more information.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for error in checking the E-M algorithm's convergence. Default: 1e-04.
#' @param zeta The initial value for the proportion of units that are not clustered. Default: 0.1. If 0, no singleton cluster is fitted.
#' @param out The file directory for writing fitting information in each E-M iteration. Default: NULL (no information is outputted).
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
MBASIC <- function(Y, S, fac, J=NULL, maxitr = 100, struct = NULL, para = NULL,  family="lognormal", method = "MBASIC", zeta = 0.1, tol = 1e-4, out = NULL, X = NULL) {

  write.out(out, "Started")
  if(! method %in% c("SE-HC", "SE-MC", "PE-MC", "MBASIC")) {
    message("Error: method must be one of SE-HC, SE-MC, PE-MC or MBASIC.")
    return
  }
  
  ## prespecified
  K <- length(unique(fac))
  I <- ncol(Y)
  N <- nrow(Y)
  if(length(fac) != N)
    message("Error: total number of replicates do not match with the number of rows in Y")
  
  ## design matrix D is K by N
  D <- matrix(0, nrow = K, ncol = N)
  for(k in 1:K) {
    D[k, fac == unique(fac)[k]] <- 1
  }
  SampleToExp <- apply(D, 2, function(x) which(x==1))
  
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
  
  numpar <- 2 * S * N + I * (S - 1) + J + sum(apply(struct, 2, function(x) length(unique(x)))) * (S - 1)
  
  outitr <- 0
  totallik <- oldlik <- 0
  alllik <- allerr <- allzeta <- bestW <- allmisclass <- matchId1 <- W.err <- matchId2 <- allari <- numeric(0)
  maxlik <- -Inf
  
  write.out(out, "Initialized parameters")
  
  ## initialize mu and variance
  Sigma <- Mu <- matrix(0, nrow = N, ncol = S)

  InitMuSigma()
  
  b <- rep(0, I)
  B <- matrix(rep(b, each = K), nrow = K)
  
  ## initialize the matrices by hierarchical clustering
  ## in constructing Z, cluster all locis
  ## This gives deterministic initialization
  ProbMat <- matrix(0, nrow = K * S, ncol = I)
  InitProbMat()

  if(method != "MBASIC") {
      ## SE-HC, PE-MC or SE-MC method
      Pi <- matrix(1/S, nrow = K, ncol = S)
      ## em step to estimate Theta
      allpar <- c(c(Mu), c(Sigma), c(Pi))
      oldpar <- 0
      for(itr in 1:maxitr) {
          ## check for convergence
          if(max(abs(oldpar - allpar))< tol)
              break
          
          ## M step
          for(s in seq_len(S)) {
              idx <- seq_len(K) + (s-1) * K
              Pi[,s] <- apply(ProbMat[idx,], 1 ,mean)
          }
          
          UpdateMuSigma()
          
          ## E step
          UpdateProbMat()
      }## finish iteration
      
      Theta <- matrix(-1, nrow = K, ncol = I)
      for(k in seq_len(K)) {
          idx <- k + K * (seq_len(S) - 1)
          Theta[k,] <- apply(ProbMat[idx,], 2, which.max)
      }
      if(!is.null(para))
          allerr <- mean(Theta != para$Theta)
      
      if(method != "PE-MC") {
          ret <- MBASIC.state(Theta, J=J, zeta = zeta, struct = struct, method = method, maxitr = maxitr, tol = tol, para = para, out = out)
          
          conv <- FALSE
          if(ret@converged & itr < maxitr)
              conv <- TRUE
          
          ## Pi is the proportion for components in the k experiment to have state s
          ## Pi is different from Z. Z is the posterior probability.
          
          return(new("MBASICFit",
                     Theta = Theta,
                     W = ret@W,
                     Z = ret@Z,
                     b = ret@b,
                     lik = ret@lik,
                     alllik = ret@alllik,
                     zeta = ret@zeta,
                     Mu = Mu,
                     Sigma = Sigma,
                     probz = ret@probz,
                     P = ret@P,
                     converged = conv,
                     Theta.err = tail(allerr, 1),
                     ARI = ret@ARI,
                     W.err = ret@W.err,
                     MisClassRate = ret@MisClassRate
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
  predZ <- Z
  b.prob <- b
  clustOrder <- .orderCluster(W, struct)
  W <- W[, clustOrder]
  Z <- Z[, clustOrder]
  W <- .structure(W, struct)
  ## initialize p, probz
  P <- matrix(0, nrow = I, ncol = S)
  if(family != "binom") {
      for(s in seq_len(S)) {
          idx <- seq_len(K) + K * (s - 1)
          P[, s] <- apply(ProbMat[idx, ], 2, mean)
      }
  } else {
      P[, 2] <- 1
  }
  probz <- apply(rbind(Z, diag(rep(1, J))), 2, mean)
  
  ## EM algorithm
  oldpar <- 0
  newpar <- c(c(W), probz, zeta, c(P), c(Mu), c(Sigma))
  
  for(outitr in seq_len(maxitr)) {
    
    if(max(abs(newpar - oldpar)) < tol) {
      break
    }
    
    ## transform everything into matrices
    B <- matrix(rep(b, each = K), nrow = K)
    
    W.lik <- matrix(0, nrow = K, ncol = J * S)
    for(s in seq_len(S)) {
        W.lik[, seq_len(J) + J * (s - 1)] <- W[seq_len(K) + K * (s - 1) ,]
    }
    
    if(outitr == 1 | method == "MBASIC") {
        ## only compute the PDF once if the method is PE-MC
        PDF <- matrix(0, nrow = N, ncol = I * S)
        for(s in seq_len(S)) {
            PDF[, seq_len(I) + I * (s - 1)] <- logdensity(Y, Mu[, s], Sigma[,s], X, family)
        }
        PDF[PDF > 5] <- 5
        PDF[PDF < -5000] <- -5000
        PDF[is.na(PDF)] <- mean(PDF, na.rm = TRUE)
        PDF <- crossprod(t(D), PDF)
    }
    
    oldlik <- totallik
    totallik <- .Call("loglik", W.lik, P, zeta, probz, PDF, package="MBASIC")
    write.out(out, paste("itr", outitr, "lik", round(totallik, 2), "zeta", round(zeta, 2)))
    ##write.out(out, paste("loglik in C =", totallik))
    
    alllik <- c(alllik, totallik)
    allzeta <- c(allzeta, zeta)
    Theta <- matrix(-1, nrow = K, ncol = I)
    for(k in seq_len(K)) {
        idx <- k + K * (seq_len(S) - 1)
        Theta[k,] <- apply(ProbMat[idx,], 2, which.max)
    }
    if(length(names(para)) != 0) {
        allerr <- c(allerr, mean(para$Theta != Theta))
        
        ## compute misclassification rate
        W.f <- matrix(0, nrow = K * S, ncol = J)
        for(s in seq_len(S))
            W.f[s + S * seq(0, K - 1),] <- W[seq_len(K) + K * (s - 1),]
        
        mc <- matchCluster(W.f, para$W, predZ, para$Z, b.prob, para$non.id)
        
        write.out(out, paste("mis-class rate ", mc$mcr))
        write.out(out, paste("Error for W ",  round(mc$W.err, 3)))
        allmisclass <- c(allmisclass, mc$mcr)
        W.err <- c(W.err, mc$W.err)
        allari <- c(allari, mc$ari)
        write.out(out, paste("ARI ", mc$ari))
        write.out(out, paste("loglik", totallik, "err", round(allerr[length(allerr)], 2)))
    }
    
    if(maxlik < totallik) {
        maxlik <- totallik
        bestb <- b.prob
        bestW <- W
        bestTheta <- Theta
    }
    
    ## E step
    ## M step for some parameters
    if(zeta > 0) {
        estep.result <- .Call("e_step", W.lik, P, zeta, probz, PDF, package = "MBASIC")
    } else {
        estep.result <- .Call("e_step1", W.lik, P, probz, PDF, package = "MBASIC")
    }
    
    ## Expected Theta matrix
    ProbMat <- estep.result[["Theta_mean"]]
    ## Maximizers
    zeta <- estep.result[["zeta"]]
    if(family != "binom") {
        P <- estep.result[["P"]]
    }
    W <- estep.result[["W"]]
    probz <- estep.result[["probz"]]
    predZ <- estep.result[["predZ"]]
    b.prob <- estep.result[["b_prob"]]
    
    W.aug <- matrix(0, nrow = K * S, ncol = J)
    for(s in 1:S) {
        W.aug[seq_len(K) + K * (s - 1),] <- W[,  seq_len(J) + J * (s - 1)]
    }
    clustOrder <- .orderCluster(W.aug, struct)
    W <- W.aug[, clustOrder]
    W <- .structure(W, struct)
    probz <- probz[clustOrder]
    predZ <- predZ[, clustOrder]
    
    ## format ProbMat
    ProbMat.Format <- matrix(0, nrow = K * S, ncol = I)
    for(s in 1:S) {
        ProbMat.Format[(s - 1) * K + seq_len(K) ,] <- ProbMat[, seq_len(I) + I * (s - 1)]
    }
    ProbMat <- ProbMat.Format
    
    if(method != "PE-MC") {
        ## skip this step if the method is PE-MC
        ## M-step for Mu and Sigma
        UpdateMuSigma()
    }
    
    ## convert everything to matrices
    B <- matrix(rep(b, each = K), nrow = K)
    
    oldpar <- newpar
    newpar <- c(c(W), probz, zeta, c(P), c(Mu), c(Sigma))
    
  }## finish outer loop
  
  conv <- FALSE
  if(outitr < maxitr)
    conv <- TRUE
  
  new("MBASICFit",
      Theta = bestTheta,
      W = bestW,
      Z = predZ,
      b = bestb,
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
      Theta.err = tail(allerr, 1),
      ARI = tail(allari, 1),
      W.err = tail(W.err, 1),
      MisClassRate = tail(allmisclass, 1),
      Struct = struct
    )
}

InitProbMat <- function() {
    Inherit()
    totalF <- matrix(0, nrow = K, ncol = I)
    F1 <- matrix(0, nrow = K * S, ncol = I)
    for(s in 1:S) {
        idx <- (s-1) * K + seq_len(K)
        F1[idx,] <-  exp(crossprod(t(D) , logdensity(Y, Mu[,s], Sigma[,s], X, family)))
        totalF <- totalF + F1[idx,]
    }
    totalF <- t(matrix(rep(c(t(totalF)), S), nrow = I))
    ProbMat <- F1 / totalF
    ProbMat[totalF == 0] <- 1/S
    maxProb <- max(na.omit(ProbMat[ProbMat != 1]))
    minProb <- min(na.omit(ProbMat[ProbMat != 0]))
    ProbMat[ProbMat > maxProb] <- max(c(0.999, na.omit(maxProb)))
    ProbMat[ProbMat < minProb] <- min(c(0.001, na.omit(minProb)))
    ProbMat[is.na(ProbMat)] <- mean(ProbMat, na.rm = TRUE)
    assign("ProbMat", ProbMat, envir = parent.frame())
}

InitMuSigma <- function() {
    Inherit()
    for(s in 1:S) {
        if(family == "lognormal") {
            Y.sec <- c(Y)[c(Y) <= quantile(c(Y), s / S) & c(Y) >= quantile(c(Y), (s - 1) / S)]
            Y.sec <- log(Y.sec + 1)
        } else if(family == "negbin") {
            Y.sec <- c(Y)[c(Y) < quantile(c(Y), s / S) & c(Y) > quantile(c(Y), (s - 1) / S)]
        } else {
            ## gamma-binomial distribution
            ratio <- Y / X
            ratio[X == 0] <- mean(Y[X > 0] / X[X > 0])
            Y.sec <- c(ratio)[ratio <= quantile(ratio, s / S) & ratio >= quantile(ratio, (s - 1) / S)]
        }
        m1 <- mean(Y.sec)
        m2 <- mean(Y.sec * Y.sec)
        MomentEstimate()
    }
    assign("Mu", Mu, envir = parent.frame())
    assign("Sigma", Sigma, envir = parent.frame())
}

UpdateMuSigma <- function() {
    Inherit()
    for(s in seq_len(S)) {
        idx <- SampleToExp + (s - 1) * K
        if(family == "lognormal") {
            m1 <- apply(log(Y + 1) * ProbMat[idx, ], 1, sum) / apply(ProbMat[idx, ], 1, sum)
            m2 <- apply(log(Y + 1) ^ 2 * ProbMat[idx, ], 1, sum) / apply(ProbMat[idx, ], 1, sum)
        } else if(family == "negbin"){
            ## negative binomial family
            m1 <- apply(Y * ProbMat[idx, ], 1, sum) / apply(ProbMat[idx, ], 1, sum)
            m2 <- apply(Y ^ 2 * ProbMat[idx, ], 1, sum) / apply(ProbMat[idx, ], 1, sum)
        } else {
            ## gamma-binomial distribution
            ratio <- Y / X
            ratio[X == 0] <- mean(Y[X > 0] / X[X > 0])
            m1 <- apply(ratio * ProbMat[idx, ], 1, sum) / apply(ProbMat[idx, ], 1, sum)
            m2 <- apply(ratio * ratio * ProbMat[idx, ], 1, sum) / apply(ProbMat[idx, ], 1, sum)
        }
        MomentEstimate()
    }
    ## order the means
    od <-  apply(Mu, 1, order)
    Mu <- matrix(Mu[cbind(rep(seq_len(N), each = S), c(od))], ncol = S, byrow = TRUE)
    Sigma <- matrix(Sigma[cbind(rep(seq_len(N), each = S), c(od))], ncol = S, byrow = TRUE)

    assign("Mu", Mu, envir = parent.frame())
    assign("Sigma", Sigma, envir = parent.frame())
}

UpdateProbMat <- function() {
    Inherit()
    F1  <- matrix(0, nrow = K * S, ncol = I)
    for(s in 1:S) {
        idx <- (s-1) * K + seq_len(K)
        F1[idx,] <-  crossprod(t(D), logdensity(Y, Mu[, s], Sigma[, s], X, family)) + log(Pi[,s])
    }
    F1[is.na(F1)] <- -5000
    F1[F1 <  -5000] <- -5000
    F.max <- t(matrix(apply(matrix(t(F1), ncol = S), 1, max), nrow = I))
    totalF <- matrix(0, nrow = K, ncol = I)
    for(s in seq_len(S)) {
        idx <- seq_len(K) + (s-1) * K
        F1[idx,] <- exp(F1[idx,] - F.max)
        totalF <- totalF + F1[idx,]
    }
    totalF <- t(matrix(rep(c(t(totalF)), S), nrow = I))
    
    ProbMat <- F1 / totalF
    maxProb <- max(ProbMat[ProbMat != 1])
    minProb <- min(ProbMat[ProbMat != 0])
    ProbMat[ProbMat > maxProb] <- maxProb
    ProbMat[ProbMat < minProb] <- minProb

    assign("ProbMat", ProbMat, envir = parent.frame())
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
  } else if(family == "negbin") {
    return(dnbinom(y, mu = mu, size = sigma, log = TRUE))
  } else {
      a <- mu / (1 - mu) * sigma
      b <- sigma
      return(
          log(beta(a + y, x - y + b)) -
              log(beta(a, b)) + log(choose(x, y))
      )
  }
}

MomentEstimate <- function() {
    for(v in c("m1", "m2", "Mu", "Sigma", "family", "s")) {
        assign(v, parent.frame()[[v]])
    }
    if(family == "lognormal") {
        Mu[, s] <- m1
        m2 <- m2 - m1 * m1
        m2[m2 < 0.01] <- 0.01
        Sigma[, s] <- m2
    } else if(family == "negbin") {
        Mu[, s] <- m1
        m2 <- m2 - m1 * m1
        m2 <- m1 / (m2 / m1 - 1)
        m2[m2 < 0] <- 100
        Sigma[, s] <- m2
    } else {
        m2 <- m2 - m1 * m1
        a.plus.b <- 1 / (m1 * (1 - m1)) - 1
        a <- m1 * a.plus.b
        b <- a.plus.b - a
        a[a < 0.01] <- 0.01
        b[b < 0.01] <- 0.01
        Mu[, s] <- a / (a + b)
        Sigma[, s] <- b
    }
    assign("Mu", Mu, envir = parent.frame())
    assign("Sigma", Sigma, envir = parent.frame())
}
