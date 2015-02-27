structure <- function(W, struct) {
  S <- ncol(W) / ncol(struct)
  K <- nrow(struct)
  oldW <- W
  for(j in 1:ncol(struct))
    for(l in unique(struct[,j])) {
      total <- 0
      for(s in 1:S) {
        idx = j + ncol(struct) * (s - 1)
        total <- total + sum(W[struct[, j] == l, idx])
      }
      if(total > 0)
        for(s in 1:S) {
          idx = j + ncol(struct) * (s - 1)
          W[struct[, j] == l, idx] <-
            sum(oldW[struct[, j] == l, idx]) / total
        }
      else
        for(s in 1:S) {
          idx = j + ncol(struct) * (s - 1)
          W[struct[, j] == l, idx] <- 1 / S
        }
    }
  return(W)
}

.structure <- function(W, struct) {
  S <- nrow(W) / nrow(struct)
  K <- nrow(struct)
  J <- ncol(W)
  oldW <- W
  for(j in seq_len(J)) {
    for(l in unique(struct[,j])) {
      total <- 0
      for(s in 1:S) {
        idx = which(struct[,j] == l) + (s - 1) * K
        W[idx, j] <- mean(W[idx, j])
      }
    }
    W[,j] <- W[, j] / rep(apply(matrix(W[,j], ncol = S), 1, sum), S)
  }
  return(W)
}

orderCluster <- function(W, struct) {
  ret <- NULL
  J <- ncol(struct)
  K <- nrow(struct)
  workSet <- seq_len(J)
  S <- ncol(W) / J
  for(j in seq_len(ncol(struct) - 1)) {
    Wtemp <- W[, workSet]
    for(s in seq_len(S - 1))
      Wtemp <- rbind(Wtemp, W[, workSet + J * s])
    if(var(struct[,j]) == 0)
      r2 <- apply(Wtemp,
                  2,
                  var
                )
    else
      r2 <- apply(Wtemp,
                  2,
                  function(x)
                  summary(lm(x~as.factor(paste(rep(1:S, each = K), rep(struct[,j] , S), sep = "|"))))$sigma
                )
    id <- workSet[which.min(r2)]
    ret <- c(ret, id)
    workSet <- setdiff(workSet, id)
  }
  ret <- c(ret, workSet)
  return(ret)
}

.orderCluster <- function(W, struct) {
  J <- ncol(struct)
  workSet <- seq_len(J)
  K <- nrow(struct)
  S <- nrow(W) / nrow(struct)
  ret <- NULL
  for(j in seq_len(ncol(struct) - 1)) {
    Wtemp <- W[, workSet]
    r2 <- apply(Wtemp,
                2,
                function(x)
                summary(lm(x~as.factor(paste(rep(1:S, each = K), rep(struct[,j] , S), sep = "|"))))$sigma
                )
    if(length(which.min(r2))==0)
      id <- workSet[1]
    else
      id <- workSet[which.min(r2)]
    ret <- c(ret, id)
    workSet <- setdiff(workSet, id)
  }
  ret <- c(ret, workSet)
  return(ret)
}

#' @name write.out
#' @title Write the message to an output file.
#' @description Write the message to an output buffer.
#' @param out Either the file name, or NULL. If NULL, the message is written to the standard output.
#' @param msg The string for the message to be written.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @export
write.out <- function(out, msg) {
    if(!is.null(out)) {
  	if(is.na(out))
	   return(NULL)
    write(msg, file = out, append=TRUE, sep = "\n")
    }  else
    message(msg)
}

logit <- function(x) {
    res <- x
    res[res>=1] <- Inf
    res[res<=0] <- -Inf
    res[res<1 & res> 0] <- log(res[res<1 & res> 0]/(1 - res[res<1 & res> 0]))
    return(res)
}

invlogit <- function(x) {
    if(x > 0)
        return(1 / (1 + exp(-x)))
    else
        return(exp(x) / (1 + exp(x)))
}

## pen <- function(W, Z, b) {    lambda1 * sum(abs(ProbMat - tcrossprod(W, Z)) * rep(1 - b, each = K)) + lambda2 * sum(b) }

group2mat <- function(x) {
        mat <- matrix(0, nrow = nrow(x), ncol = nrow(x))
        for(i in 1:ncol(x)) {
                id <- which(x[, i] == 1)
                mat[id, id] <-  1
        }
        return(mat)
}

matchCluster <- function(W, W.true, predZ, Z.true, b.prob, non.id) {
  if(is.null(b.prob))
    b.prob <- rep(0, nrow(predZ))
  
  distW <- matrix(0, nrow = ncol(W), ncol = ncol(W.true))
  for(j1 in seq_len(ncol(W))) 
    for(j2 in seq_len(ncol(W.true)))
      distW[j1, j2] <- sum((W[, j1] - W.true[, j2]) ^ 2)
  
  matchId1 <- matrix(-1, nrow = ncol(W), ncol = 2)
  matchId2 <- matrix(-1, nrow = ncol(W.true), ncol = 2)
  
  numMissClass <- 0
  Z.pred <- cbind(predZ * (1- b.prob), b.prob)
  Z.true <- cbind(Z.true, 0)
  Z.true[non.id,] <- 0
  Z.true[non.id, ncol(W.true) + 1] <- 1
  
  for(j in seq_len(ncol(W))) {
    matchId1[j, 1] <- which.min(distW[j,])
    matchId1[j, 2] <- j
    numMissClass <- numMissClass + sum(Z.pred[, matchId1[j, 2]] * (1 - Z.true[, matchId1[j, 1]]))
  }
  
  for(j in seq_len(ncol(W.true))) {
    matchId2[j, 1] <- j
    matchId2[j, 2] <- which.min(distW[, j])
    numMissClass <- numMissClass + sum((1 - Z.pred[, matchId2[j, 2]]) * Z.true[, matchId2[j, 1]])
  }
  
  numMissClass <- numMissClass + sum(Z.pred[, ncol(W) + 1] * (1 - Z.true[, ncol(W.true) + 1]))
  numMissClass <- numMissClass + sum((1 - Z.pred[, ncol(W) + 1]) * Z.true[, ncol(W.true) + 1])

  mcr <- numMissClass / 2/ nrow(Z.true)
  
  W.err <- sqrt(
                (
                 sum((W[, matchId1[,2]] -W.true[, matchId1 [,1]]) ^ 2) +
                 sum((W[, matchId2[,2]] -W.true[, matchId2 [,1]]) ^ 2)
               ) / (ncol(W) + ncol(W.true)) / nrow(W.true)
              )

  ari <- adjustedRandIndex(
                           apply(Z.pred, 1, function(x) which.max(x)),
                           apply(Z.true, 1, function(x) which.max(x))
                         )

  return(list(mcr = mcr, W.err = W.err, ari = ari, matchId1 = matchId1, matchId2 = matchId2))
                
}

extract <- function(x, head, tail) {
  m <- regexpr(head, x)
  tmp <- regmatches(x, m, invert = TRUE)
  res <- NULL
  for(i in seq_along(tmp)) {
    res <- c(res,
             regmatches(tmp[[i]][2], regexpr(tail, tmp[[i]][2]))
           )
  }
  return(res)
}

trimLogValue <- function(x) {
  x[x > 5] <- 5
  x[x < -5000] <- -5000
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}

trimProbValue <- function(x) {
  maxProb <- max(na.omit(x[x != 1]))
  minProb <- min(na.omit(x[x != 0]))
  x[x > maxProb] <- max(c(0.999, na.omit(maxProb)))
  x[x < minProb] <- min(c(0.001, na.omit(minProb)))
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}

estep <- function(W, P, V, zeta, probz, PDF, designMap, stateMap, unitMap) {
  I <- ncol(PDF)
  J <- ncol(W)
  K <- ncol(designMap)
  N <- nrow(designMap)
  S <- ncol(stateMap)
  M <- nrow(stateMap)

  PDF <- exp(trimLogValue(PDF))

  ## compute joint PDF
  tmpPDF <- PDF * matrix(c(V), nrow = N * M, ncol = I)
  tmpPDF <- crossprod(unitMap, tmpPDF) ## NS x I
  jointPDF <- matrix(0, nrow = K * S, ncol = I)
  for(s in seq(S)) {
    id1 <- seq(N) + (s - 1) * N
    id2 <- seq(K) + (s - 1) * K
    jointPDF[id2, ] <- crossprod(designMap, log(tmpPDF[id1, ]))
  }
  jointPDF <- exp(trimLogValue(jointPDF))

  ## compute FP
  tmp <- matrix(rep(c(t(P)), each = K), nrow = K * S, ncol = I)
  tmp <- jointPDF * tmp
  FP <- tmp[seq(K), ]
  for(s in seq(S)[-1]) {
    FP <- FP + tmp[seq(K) + (s - 1) * K, ]
  }

  ## FW
  workMat <- diag(rep(1, K))
  for(s in seq(S)[-1]) {
    workMat <- rbind(workMat, diag(rep(1, K)))
  }
  FW <- matrix(0, nrow = K * J, ncol = I)
  for(j in seq(J)) {
    FW[seq(K) + (j - 1) * K, ] <- crossprod(workMat, jointPDF * W[, j])
  }

  ## FV
  FV <- PDF * matrix(c(V), nrow = N * M, ncol = I)
  FV <- crossprod(unitMap, FV)

  FW <- trimLogValue(log(FW))
  FV <- trimLogValue(log(FV))
  FP <- trimLogValue(log(FP))
  jointPDF <- trimLogValue(log(jointPDF))
  PDF <- trimLogValue(log(PDF))

  workMat <- matrix(0, nrow = K * J, ncol = J)
  workMat[cbind(seq(K * J), rep(seq(J), each = K))] <- 1

  Wsum <- crossprod(workMat, FW)
  Wsum <- exp(trimLogValue(Wsum)) * probz * (1 - zeta)

  Psum <- exp(trimLogValue(apply(FP, 2, sum))) * zeta
  
  Zcond <- t(Wsum)
  Zcond <- Zcond / apply(Zcond, 1, sum)

  Znew <- t(Wsum) + Psum * rep(probz, each = I)
  Znew <- Znew / apply(Znew, 1, sum)

  bnew <- Psum / (Psum + apply(Wsum, 2, sum))

  Thetab <- jointPDF
  for(s in seq(S)) {
    Thetab[seq(K) + (s - 1) * K, ] <- Thetab[seq(K) + (s - 1) * K, ] - FP + rep(log(bnew) + log(P[, s]), each = K)
  }

  Thetazb <- matrix(0, nrow = K * S * J, ncol = I)
  for(j in seq(J)) {
    for(s in seq(S)) {
      Thetazb[seq(K) + (s - 1) * K + (j - 1) * S * K, ] <- jointPDF[seq(K) + (s - 1) * K, ] + log(W[seq(K) + (s - 1) * K, j]) - FW[seq(K) + K * ( j - 1), ] + rep(log(Zcond[, j]) + log(1 - bnew), each = K)
    }
  }

  Wnew <- matrix(apply(exp(trimLogValue(Thetazb)), 1, sum), ncol = J)
  for(k in seq(K)) {
    idx <- k + (seq(S) - 1) * K
    Wnew[idx, ] <- Wnew[idx, ] / rep(apply(Wnew[idx, ], 2, sum), each = S)
  }

  Pnew <- t(exp(trimLogValue(Thetab)))
  workMatrix <- matrix(0, nrow = S, ncol = K * S)
  workMatrix[cbind(rep(seq(S), each = K), seq(K * S))] <- 1
  Pnew <- tcrossprod(Pnew, workMatrix)
  Pnew <- Pnew / apply(Pnew, 1, sum)

  Thetanew <- exp(trimLogValue(Thetab))
  for(j in seq(J)) {
    Thetanew <- Thetanew + exp(trimLogValue(Thetazb[seq(K * S) + K * S * (j - 1), ]))
  }
  
  Thetanu <- matrix(0, nrow = M * N, ncol = I)
  tmp <- matrix(0, nrow = N * S, ncol = I)
  for(s in seq(S)) {
    tmp[seq(N) + N * (s - 1), ] <- designMap %*% Thetanew[seq(K) + K * (s - 1), ]
  }
  Nu <- Thetanu <- unitMap %*% tmp
  Thetanu <- Thetanu * c(V)
  Thetanu <- trimLogValue(log(Thetanu)) + PDF
  for(n in seq(N)) {
    id1 <- n + (seq(M) - 1) * N
    id2 <- n + (seq(S) - 1) * N
    Thetanu[id1, ] <- Thetanu[id1, ] - stateMap %*% FV[id2, ]
  }
  Thetanu <- exp(trimLogValue(Thetanu))

  Nu <- (1 - Nu) * c(V)
  Nu <- Nu + Thetanu
  
  Vnew <- matrix(apply(Nu, 1, sum), ncol = M)
  Vnew <- Vnew / (Vnew %*% stateMap %*% t(stateMap))
  
  return(list(
      P = Pnew,
      zeta = mean(bnew),
      probz = apply(Znew, 2, mean),
      W = Wnew,
      Theta = Thetanew,
      Theta_nu = Thetanu,
      V = Vnew,
      b_prob = bnew,
      Zcond = Zcond,
      Z = Znew,
  jointPDF = jointPDF,
  FW = FW,
  FV = FV,
  FP = FP)
         )
}
