
structure <- function(W, struct) {
  S <- ncol(W) / ncol(struct)
  K <- nrow(struct)
  oldW <- W
  for(j in 1:ncol(struct))
    for(l in unique(struct[,j])) {
      total <- 0
      for(s in 1:S) {
        idx = j + ncol(struct) * (s - 1)
        total <- total + sum(W[ struct[, j] == l, idx ])
      }
      if(total > 0)
        for(s in 1:S) {
          idx = j + ncol(struct) * (s - 1)
          W[ struct[, j ] == l, idx ] <-
            sum(oldW[ struct[, j ] == l, idx ]) / total
        }
      else
        for(s in 1:S) {
          idx = j + ncol(struct) * (s - 1)
          W[ struct[, j ] == l, idx ] <- 1 / S
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
        idx = which(struct[ ,j ] == l) + (s - 1) * K
        W[ idx, j ] <- mean(W[ idx, j ])
      }
    }
    W[ ,j ] <- W[ , j ] / rep(apply(matrix(W[ ,j ], ncol = S), 1, sum), S)
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
    Wtemp <- W[ , workSet ]
    for(s in seq_len(S - 1))
      Wtemp <- rbind(Wtemp, W[ , workSet + J * s ])
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
    id <- workSet[ which.min(r2) ]
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
    Wtemp <- W[ , workSet ]
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
    if(length(which.min(r2))==0)
      id <- workSet[ 1 ]
    else
      id <- workSet[ which.min(r2) ]
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

logdensity <- function(y, mu, sigma, family) {
  if(family == "lognormal") {
    y <- log(y + 1)
    return(-(y - mu) ^ 2 / sigma / 2 - log(sigma) / 2 - log(2 * pi) / 2 )
  } else if(family == "negbin") {
    return(dnbinom(y, mu = mu, size = sigma, log = TRUE))
  } else {
    return(dbinom(y, size = sigma, prob = mu, log = TRUE))
  }
}

logit <- function(x) {
    res <- x
    res[res>=1] <- Inf
    res[res<=0] <- -Inf
    res[res<1 & res> 0 ] <- log(res[res<1 & res> 0 ]/(1 - res[res<1 & res> 0 ]))
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
                id <- which(x[ , i ] == 1)
                mat[ id, id ] <-  1
        }
        return(mat)
}

matchCluster <- function(W, W.true, predZ, Z.true, b.prob, non.id) {
  if(is.null(b.prob))
    b.prob <- rep(0, nrow(predZ))
  
  J <- ncol(W)
  distW <- matrix(0, nrow = J, ncol = ncol(W.true))
  for(j1 in seq_len(J)) 
    for(j2 in seq_len(ncol(W.true)))
      distW[ j1, j2 ] <- sum((W[ , j1 ] - W.true[ , j2 ]) ^ 2)
  
  matchId1 <- matrix(-1, nrow = J, ncol = 2)
  matchId2 <- matrix(-1, nrow = ncol(W.true), ncol = 2)
  
  numMissClass <- 0
  Z.pred <- cbind(predZ * ( 1- b.prob), b.prob)
  Z.true <- cbind(Z.true, 0)
  Z.true[ non.id, ] <- 0
  Z.true[ non.id, ncol(W.true) + 1] <- 1
  
  for(j in seq_len(J)) {
    matchId1[ j, 1 ] <- which.min(distW[ j, ])
    matchId1[ j, 2 ] <- j
    numMissClass <- numMissClass + sum(Z.pred[ , matchId1[ j, 2 ] ] * (1 - Z.true[ , matchId1[ j, 1 ] ]))
  }
  
  for(j in seq_len(ncol(W.true))) {
    matchId2[ j, 1 ] <- j
    matchId2[ j, 2 ] <- which.min(distW[ , j ])
    numMissClass <- numMissClass + sum((1 - Z.pred[ , matchId2[ j, 2 ] ]) * Z.true[ , matchId2[ j, 1 ] ])
  }
  
  numMissClass <- numMissClass + sum(Z.pred[ , J + 1 ] * (1 - Z.true[ , ncol(W.true) + 1 ]))
  numMissClass <- numMissClass + sum((1 - Z.pred[ , J + 1 ]) * Z.true[ , ncol(W.true) + 1 ])

  mcr <- numMissClass / 2/ nrow(Z.true)
  
  W.err <- sqrt(
                (
                 mean((W[, matchId1[ ,2] ] -W.true[ , matchId1 [,1] ]) ^ 2) +
                 mean((W[, matchId2[ ,2] ] -W.true[ , matchId2 [,1] ]) ^ 2)
                ) / 2
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
