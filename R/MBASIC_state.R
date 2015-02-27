#' @name MBASIC.state
#' @title Bayesian clustering model for a state-space matrix.
#'
#' @description This function clusters a state-space matrix.
#' @param Theta A K by I matrix. The (k,i)-th entry is the state of the i-th unit under condition k. Notice that the sorted distinct values of entries in this matrix must be 1,2,...,S, where S is the total number of states.
#' @param struct A K by J matrix indicating the structures of each cluster.
#' @param J The number of clusters to be identified.
#' @param method A string for the fitting method, 'SE-HC' or 'SE-MC'(default).
#' @param para A list object that contains the true model parameters. Default: NULL. See details for more information.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for error in checking the E-M algorithm's convergence. Default: 1e-04.
#' @param zeta The initial value for the proportion of units that are not clustered. Default: 0.1.
#' @param out The file directory for writing fitting information in each E-M iteration. Default: NULL (no information is outputted).
#' @details
#'  The 'method' argument determines what fitting method will be used. The default is 'em', where an E-M algorithm is used for clustering. If 'naive', then hierarchical clustering is used.\cr
#' The 'para' argument takes a list object that is supposed to include the following fields:
#'\tabular{ll}{
#' W \tab A K by (J*S) matrix. The (k,J*(s-1)+j)-th entry is the probability that the units in cluster j has state s in the k-th experiment.\cr
#' Z \tab An I by J matrix. The (i,j)-th entry is the indicator whether the i-th unit belongs to cluster j.\cr
#' non.id \tab A binary vector of length I. The i-th entry is the indicator whether the i-th unit does not belong to any cluster.
#' }
#' This argument is intended to carry the true parameters in simulation studies. If it is not null, then the model also computes a number of metrics that describes the error in model fitting. Users should be cautious that the order of the rows and columns of matrices in the fields of para should match the Y matrix.
#' @return An object of class 'MBASICFit'.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' state.sim <- MBASIC.sim.state(I = 1000, K = 10, J = 4, S = 3, zeta = 0.1)
#' state.sim.fit <- MBASIC.state(Theta = state.sim$Theta, J = 4, method = "SE-MC", zeta = 0.1, maxitr = 100, tol = 1e-04)
#' @useDynLib MBASIC
#' @export
MBASIC.state <- function(Theta, J, struct = NULL, method = "SE-MC", zeta = 0.1, maxitr = 100, tol = 1e-4, para = NULL, out = NULL) {

  mc <- NULL
  if(!method %in% c("SE-HC", "SE-MC")) {
    message("method must be SE-HC or SE-MC")
    return
  }
  
  S <- max(Theta)
  if(prod(sort(unique(c(Theta))) == seq_len(S)) == 0) {
    message("The sorted distinct entries in the matrix Theta must be 1,2,...,S")
    return
  }
  
  K <- nrow(Theta)
  I <- ncol(Theta)

  ProbMat <- matrix(0, nrow = K * S, ncol = I)
  ProbMat[cbind((0:(K-1)) * S + c(Theta), rep(1:ncol(Theta), each = K))] <- 1
  
  ## use hierarchical clustering
  d <- .Call("hamming", Theta, package = "MBASIC")
  d <- d + t(d)
  
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
  
  alllik <- allerr <- allwerr <- allmisclass <- P <- NULL
  
  if(method == "SE-HC") {
    d <- as.dist(d)
    fit <- hclust(d)
    groups <- cutree(fit, k = J)
    Z <- matrix(0, nrow = I, ncol = J)
    Z[cbind(1:I, groups)] <- 1
    b.prob <- rep(0, I)
    
    W <- tcrossprod(ProbMat, t(Z))
    W <- W / rep(apply(matrix(W, nrow = S), 2, sum), each = S)

    if(!is.null(para)) {
       mc <- matchCluster(W, para$W, Z, para$Z, b.prob, para$non.id, S)
       allwerr <- mc$W.err
       allmisclass <- mc$mcr
    } 
    Zcond <- predZ <- Z
    probz <- apply(Z, 2, mean)
    
  } else {

    ## iteratively using the EM algorithm

    ## initialize groups
    
    mind <- apply(d, 1, function(x) min(x[x > 0]))
    thr <- quantile(mind, 1 - zeta)
    id <- which(mind <= thr)
    b <- rep(1, I)
    b[id] <- 0
    d <- .Call("hamming", Theta[, id], package="MBASIC")
    d <- as.dist(d + t(d))
    fit <- hclust(d)
    groups <- cutree(fit, k = J)
    
    Z <- matrix(0, nrow = I, ncol = J)
    Z[cbind(1:I, sample(1:J, I, replace = TRUE))] <- 1
    Z[id,] <- 0
    Z[cbind(id, groups)] <- 1

    W <- matrix(0, nrow = K * S, ncol = J)
    for(j in seq_len(J)) {
      if(sum(groups == j) > 0) {
        W[,j] <- apply(ProbMat[, id[groups == j]], 1, mean)
      }
    }

    W[W < tol] <- tol
    W[W > 1-tol] <- 1 - tol

    P <- matrix(0, nrow = I, ncol = S)
    for(i in seq_len(I))
      P[i,] <- apply(matrix(ProbMat[, i], nrow = S), 1, mean)

    oldW <- W
    
    W  <- matrix(0, nrow = K, ncol = J * S)
    PDF <- matrix(0, nrow = K, ncol = I * S)
    for(s in seq_len(S)) {
      W[, seq_len(J) + J * (s - 1)] <- oldW[s + S * (seq_len(K) - 1), ]
      PDF[, seq_len(I) + I * (s - 1)] <- ProbMat[s + S * (seq_len(K) - 1), ]
    }

    probz <- apply(Z, 2, mean)

    allpar <- c(c(W), c(P), zeta, probz)
    
    for(outitr in seq_len(maxitr)) {

      oldpar <- allpar
      
      totallik <- .Call("loglik_theta", W, P, zeta, probz, PDF, package = "MBASIC")
      alllik <- c(alllik, totallik)
      
      mcmc.result <- .Call("e_step_theta", W, P, zeta, probz, PDF, package = "MBASIC")
      
      ## Maximizers
      zeta <- mcmc.result[["zeta"]]
      W <- mcmc.result[["W"]]
      ##W <- W[, 1:J] / (W[, 1:J] + W[, (1:J) + J])
      probz <- mcmc.result[["probz"]]
      predZ <- mcmc.result[["Z"]]
      Zcond <- mcmc.result[["Zcond"]]
      b.prob <- mcmc.result[["b_prob"]]
      
      W.aug <- matrix(0, nrow = K * S, ncol = J)
      for(s in 1:S) {
        W.aug[seq_len(K) + K * (s - 1),] <- W[,  seq_len(J) + J * (s - 1)]
      }
      clustOrder <- .orderCluster(W.aug, struct)
      W.aug <- W.aug[, clustOrder]
      W.aug <- .structure(W.aug, struct)
      probz <- probz[clustOrder]
      predZ <- predZ[, clustOrder]
      Zcond <- Zcond[, clustOrder]
      
      if(!is.null(para)) {
        Z.format <- matrix(0, nrow = I, ncol = J)
        Z.format[cbind(1:I, apply(predZ, 1, which.max))] <- 1
        W.format <- matrix(0, nrow = K * S, ncol = J)
        for(s in 1:S) {
          idx <- s + S * seq(0, K - 1)
          W.format[idx,] <- W.aug[seq_len(K) + K * (s - 1),]
        }
        mc <- matchCluster(W.format, para$W, Z.format, para$Z, b.prob, para$non.id, S)
        allwerr <- c(allwerr, mc$W.err)
        allmisclass <- c(allmisclass, mc$mcr)
      }

      allpar <- c(c(W), c(P), zeta, probz)

      if(max(allpar - oldpar) < tol)
        break
      
    }
    
  }

  conv <- TRUE
  if(method == "SE-MC") {
    if(outitr >= maxitr)
      conv <- FALSE
  }
  
  mcr <- werr <- ari <- numeric(0)
  if(!is.null(mc)) {
    mcr <- mc$mcr
    werr <- mc$W.err
    ari <- mc$ari
  }

  write.out(out, "finished fitting Theta")

  if(method == "SE-MC") {
    return(new("MBASICFit",
               W = W,
               Z = predZ,
               clustProb = cbind(b.prob, Zcond * (1 - b.prob)),
               b = b.prob,
               lik = tail(alllik, 1),
               alllik = alllik,
               zeta = zeta,
               probz = probz,
               P=P,
               converged = conv,
               MisClassRate = mcr,
               W.err = werr,
               ARI = ari,
               Struct = struct
               ))
  } else {
    return(new("MBASICFit",
               W = W,
               Z = predZ,
               clustProb = cbind(b.prob, Zcond * (1 - b.prob)),
               b = b.prob,
               zeta = 0,
               probz = probz,
               converged = conv,
               MisClassRate = mcr,
               W.err = werr,
               ARI = ari,
               Struct = struct
               ))
  }
  
}
