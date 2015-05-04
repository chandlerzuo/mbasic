#' @name MBASIC.binary
#' @title Bayesian clustering model for binary state matrix with prior estimated background means. This function is deprecated.
#' @param Y An N by I matrix containing the data from N experiments across I observation units (loci).
#' @param Mu0 An N by I matrix for the prior estimated mean for the background state, for N experiments across the I observation units (loci).
#' @param fac A vector of length N denoting the experimental condition for each replicate.
#' @param struct A matrix indicating the levels of the signal matrix.
#' @param J The number of clusters to be identified.
#' @param family The distribution of family to be used. Either "lognormal" or "negbin". See details for more information.
#' @param burnin An integer value for the number of iterations in initialization. Default: 20.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for error in checking the E-M algorithm's convergence. Default: 1e-04.
#' @param nsig The number of mixture components for the distribution of the signal state.
#' @param zeta The initial value of the proportion of unclustered units. Default: 0.2.
#' @param out The file directory for writing fitting information in each E-M iteration. Default: NULL (no information is outputted).
#' @param outfile The file directory for writing the intermediate results every 10 E-M iterations. This can be useful when the running time until final convergence is long. Default: NULL (no intermediate result is saved).
#' @param init.mod A 'MBASICFit' class object taken as the initial model to initialize parameters. This object can be an intermediate result from a not yet converged E-M algorithm, or a fitted model with smaller number of clusters. The user must be cautious in providing this initial model, since it must be fitted using the same data as the argument of this function.
#' @details
#' Function MBASIC.binary currently supports two different distributional families: log-normal and negative binomial. This should be specified by the 'family' argument.\cr
#' For the log-normal distributions, log(Y+1) is modeled as normal distributions. For experiment n, if locus i is unenriched, distribution for log(Y[n,i]+1) is N(e[n] * Mu0[n,i], sigma0[n]). If locus i is enriched and the enrichment state is s, the distribution of log(Y[n,i]+1) is N(mu1[n,s], sigma1[n,s]).\cr
#' For the negative binomial distributions, the meanings of mu1, sigma1, sigma0 changes. For experiment n, if locus i is unenriched, distribution of Y[n,i] is NB(Mu0[n,i] * e[n], sigma0[n]). Otherwise, if locus i is enriched and the enrichment component is s, the distribution of Y[n]-1 is NB(mu1[n,s], sigma[n,s]). In this package, NB(mu, a) denotes the negative-binomial distribution with mean mu and size a (i.e. the variance is mu*(1+mu/a)).
#' @useDynLib MBASIC
#' @return A 'MBASICFit' class object.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' \dontrun{
#' ## Simulate a data
#' dat <- MBASIC.sim.binary(I = 100, fac = rep(1:5, each = 2), J = 3, f = 5, family = "lognormal")
#' ## Fit the model
#' fit <- MBASIC.binary(t(dat$Y), Mu0 = t(dat$Mu0), fac = rep(1:5, each = 2), J=3, struct = NULL, family="lognormal")
#' }
#' @export
MBASIC.binary <- function(Y, Mu0, fac, J=NULL, zeta=0.2, maxitr = 100, burnin = 20, outfile=NULL, out=NULL, init.mod = NULL, struct = NULL, family="lognormal", tol = 1e-4, nsig = 1, min.count = 5) {
  
    write.out(out, "starting...")
    ## prespecified
    K <- length(unique(fac)) 
    I <- ncol(Y)
    N <- nrow(Y)
    if(length(fac) != N)
      stop("Error: total number of replicates do not match with the number of rows in Y")

    if(prod(dim(Y) == dim(Mu0)) != 1)
      stop("Error: dimensions for Y and Mu0 must be the same.")
    
    ## index
    ##id <- rep(n[1], K)
    ##for(k in 2:K) {
    ##  id[k] <- id[k-1] + n[k]
    ##}

    if(!is.null(init.mod) & class(init.mod) != "MBASICFit")
      stop("Error: init.mod must be an object of class 'MBASICFit'.")

    ## design matrix D is K by N
    D <- matrix(0, nrow = K, ncol = length(fac))
    for(k in 1:K) {
      D[ k, fac == unique(fac)[ k ] ] <- 1
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

    numpar <- (2 + 2 * nsig) * N + (nsig - 1) * K + 1 + I + (J - 1) + sum(apply(struct, 2, function(x) length(unique(x))))

    outitr <- 0
    totallik <- oldlik <- 0
    alllik <- bestW <- bestV<- bestzeta <- NULL
    maxlik <- -Inf

    lambda <- logit(1 - zeta)

    ## normalize the Mu0
    if(family == "lognormal") {
      Mu0 <- Mu0 * rep(apply(log(Y + 1), 1, mean) / apply(Mu0, 1, mean), ncol(Mu0))
    } else {
      Mu0 <- Mu0 * rep(apply(Y, 1, mean) / apply(Mu0, 1, mean), ncol(Mu0))
    }
    
    if(is.null(init.mod)) {
      ## initialize
      e <- rep(0.8, N)
      prior.bind <- 0.5

      if(family == "lognormal") {
        sigma1 <- sigma0 <- rep(0.5, N)
        mu1 <- quantile(c(log(Y + 1) + 1), (2 - prior.bind) / 2) / quantile(c(log(Y + 1) + 1), (1 - prior.bind) / 2) * apply(Mu0, 1, mean) * e
      } else {
        v <- var(c(Y)[ Y > quantile(c(Y), 1 - prior.bind)])
        mu1 <- quantile(c(Y+1), (2 - prior.bind) / 2) / c(quantile(c(Y+1), (1 - prior.bind) / 2)) * apply(Mu0, 1, mean) * e
        m1 <- mean(mu1)
        if(v > mean(m1)) {
          sigma1 <- rep(m1 * m1 / (v - m1), N)
        }      else{
          sigma1 <- rep(100, N)
        }
        v <- var(c(Y)[ Y < median(c(Y))])
        if(is.na(v)) v <- 1
        m1 <- mean(Mu0) * mean(e)
        if(v > m1) {
          sigma0 <- rep(m1 * m1 / (v - m1), N)
        }  else{
          sigma0 <- rep(100, N)
        }
        mu1[ mu1 < 1 ] <- 1
      }
      mu1 <- matrix(mu1, ncol = nsig, nrow = N) * rep(seq_len(nsig), each = N)
      sigma1 <- matrix(sigma1, ncol = nsig, nrow = N)

      ## initialize the distribution for each replicate
      write.out(out, "Initialize the distribution for each replicate")
      p <- cbind(1 - prior.bind,matrix(prior.bind / nsig, nrow = K, ncol = nsig))
      allpar <- c(c(mu1), c(sigma1), e, sigma0, c(p))
      oldpar <- 0
      for(inititr in seq_len(burnin)) {
        if(max(abs(allpar - oldpar)) < tol)
          break

        Sigma0 <- matrix(rep(sigma0, I), nrow = N)
        E <- matrix(rep(e, I), nrow = N)
        totalF <- matrix(0, nrow = K, ncol = I)
        F1 <- matrix(0, nrow = K ,  ncol = I * (nsig + 1))
        if(family == "lognormal") {
          tmpF <- -(log(Y + 1) - Mu0 * e) ^ 2 / 2 / Sigma0
        } else {
          tmpF <- log(dnbinom(Y , mu = Mu0 * e , size = Sigma0))
        }
        tmpF[ tmpF == -Inf ] <- min(tmpF[ tmpF > -Inf ], na.rm = TRUE)
        tmpF[ tmpF >=5 ] <- 5
        tmpF[ is.na(tmpF) ] <- mean(na.omit(tmpF))
        F1[ , seq_len(I) ] <- exp(crossprod(t(D), tmpF)) * p[ , 1 ]
        totalF <- totalF + F1[ , seq_len(I) ]

        for(s in seq_len(nsig)) {
          idx <- s * I + seq_len(I)
          Mu1 <- matrix(rep(mu1[ , s ], I), nrow = N)
          Sigma1 <- matrix(rep(sigma1[ , s ], I), nrow = N)
          if(family == "lognormal") {
            tmpF <- -(log(Y + 1) - Mu1) ^ 2 / 2 / Sigma1
          } else {
            tmpF <- log(dnbinom(Y - min.count , mu = Mu1 - min.count , size = Sigma1))
          }
          tmpF[ tmpF == -Inf ] <- min(tmpF[ tmpF > -Inf ], na.rm = TRUE)
          tmpF[ tmpF >= 5 ] <- 5
          tmpF[ is.na(tmpF) ] <- mean(na.omit(tmpF))
          F1[ , idx ] <- exp(crossprod(t(D), tmpF) 
                           ) * p[ , s + 1 ]
          totalF <- totalF + F1[ , idx ]
        }
        ## rbind totalF S times, creating an (K * S) * I matrix
        totalF <- matrix(rep(c(totalF), nsig + 1), nrow = K)
        ProbMat <- F1 / totalF
        ProbMat[ totalF == 0 ] <- 1 / (nsig + 1)
        maxProb <- max(ProbMat[ ProbMat != 1 ])
        minProb <- min(ProbMat[ ProbMat != 0 ])
        ProbMat[ ProbMat > maxProb ] <- max(c(0.999, maxProb))
        ProbMat[ ProbMat < minProb ] <- min(c(0.001, minProb))

        ## M step

        if(family == "lognormal") {
            ProbTheta <- crossprod(D, ProbMat[ , seq_len(I) ])
            e <- apply(log(Y + 1) * ProbTheta, 1, sum) / apply(Mu0 * ProbTheta, 1, sum)
            e[ e > 1 ] <- 1
            e[ e < 0.1 ] <- 0.1
            e[ e<0.1 ] <- 0.1
            sigma0 <- apply((log(Y + 1) - Mu0 * e) ^ 2 * ProbTheta, 1, sum) / apply(ProbTheta, 1, sum)
            for(s in seq_len(nsig)) {
                ProbTheta <- crossprod(D, ProbMat[ , s * I + seq_len(I) ])
                m0 <- apply(ProbTheta, 1, sum)
                m1 <- apply(log(Y + 1) * ProbTheta, 1, sum) / m0
                m2 <- apply(log(Y + 1) * log(Y + 1) * ProbTheta, 1, sum) / m0
                mu1[ , s ] <- m1
                sigma1[ , s ] <- m2 - m1 * m1
                sigma1[ sigma1[ , s ] <= 0, s ] <- min(sigma1[ sigma1[ , s ] > 0, s ])
            }
        } else {
            ProbTheta <- crossprod(D, ProbMat[ , seq_len(I) ])
            e <- apply(Y * ProbTheta, 1, sum) / apply(Mu0 * ProbTheta, 1, sum)
            e[ e>1 ] <- 1
            e[ e < 0.1 ] <- 0.1
            m2 <- apply(Y * Y * ProbTheta, 1, sum) / apply(ProbTheta, 1, sum)
            m1 <-  apply(Mu0 * ProbTheta * e, 1, mean)
            sigma0 <- m1 * m1 / (m2 - m1 * m1 - m1)
            sigma0[ sigma0 <= 0 ] <- 100
            for(s in seq_len(nsig)) {
                ProbTheta <- crossprod(D, ProbMat[ , s * I + seq_len(I) ])
                m0 <- apply(ProbTheta, 1, sum)
                m1 <- apply((Y - min.count) * ProbTheta, 1, sum) / m0
                m2 <- apply((Y - min.count) * (Y - min.count) * ProbTheta, 1, sum) / m0
                m1[ m1 < 0 ] <- 0
                mu1[ , s ] <- m1 + min.count
                sigma1[ , s ] <- m1 * m1 / (m2 - m1 * m1 - m1)
                sigma1[ sigma1[ , s] <= 0, s ] <- 100
            }
        }
        
        p <- t(matrix(apply(matrix(t(ProbMat), nrow = I), 2, sum), nrow = nsig + 1))
        p <- p / apply(p, 1, sum)

        lower.b <- apply(Mu0 * e, 1, mean)
        if(family == "negbin") {
            lower.b <- lower.b + 1
        }
        if(family != "negbin")## changed for using empirical Mu0
          mu1 <- (mu1 < lower.b) * lower.b + (mu1 >= lower.b) * mu1
        for(s in seq_len(nsig)[ -1 ])
          mu1[ mu1[ , s ] <= mu1[ , s-1 ], s ] <- mu1[ mu1[ , s ] <= mu1[ , s-1], s-1 ] + 1

        oldpar <- allpar
        allpar <- c(c(mu1), c(sigma1), e, sigma0, c(p))

      }

      b <- rep(0, I)
      B <- matrix(rep(b, each = K), nrow = K)
      ## initialize the matrices by hierarchical clustering
      ## in constructing Z, cluster all locis
      ## This gives deterministic initialization

      ## initialize W, Z, b
      ##if(!is.null(outfile)) {
        ##save(list = ls(), file = outfile)
      ##}
      if(FALSE) {
        write.out(out, "Initialize the distance matrix.")
        distMat <- dist(t(ProbMat[ , seq_len(I) ]))
        mind <- apply(as.matrix(distMat), 1, function(x) min(x[x>0]))
        thr <- quantile(mind, invlogit(lambda))
        id <- which(mind < thr)
      }
      unitVar <- apply(ProbMat[, seq_len(I)], 2, var)
      unitVarRank <- rank(-unitVar) / I
      id <- which(unitVarRank < zeta)
      b <- rep(1, I)
      b[id] <- 0
      if(FALSE) {
        write.out(out, "Use hierarchical clustering to generate clusters.")
        distMat <- dist(t(ProbMat[ seq_len(K), id ]))
        clustFit <- hclust(distMat)
        groups <- cutree(clustFit, k = J)
      }
      write.out(out, "Use K-means to initialize the clusters.")
      ## K-means
      clustFit <- kmeans(t(ProbMat[seq_len(K), id]), centers = J)
      groups <- clustFit$cluster
      
      Z <- matrix(0, nrow = I, ncol = J)
      Z[ cbind(1:I, sample(1:J, I, replace = TRUE))] <- 1
      Z[ id, ] <- 0
      Z[ cbind(id, groups) ] <- 1
      W <- matrix(0, nrow = K, ncol = J)
      for(j in seq_len(J)) {
        if(sum(groups == j) > 0) {
          for(k in seq_len(K)) {
            W[ k, j ] <- mean(ProbMat[ k, id[ groups == j ] ])
          }
        }
      }
      predZ <- Z
      b.prob <- b
      clustOrder <- .orderCluster(W, struct)
      W <- W[ , clustOrder ]
      Z <- Z[ , clustOrder ]
      W <- .structure(rbind(W, 1 - W), struct)[ 1:K, ]
      ## initialize p, probz
      ## p is I 
      p <- matrix(apply(ProbMat, 2, mean), nrow = I)[ , 1 ]
      probz <- apply(rbind(Z, diag(rep(1, J))), 2, mean)
      V <- t(matrix(apply(matrix(t(ProbMat[ , -seq_len(I) ]), nrow = I), 2, sum), nrow = nsig))
      V <- V / apply(V, 1, sum) ## (K * nsig)
  } else {
      mu1 <- init.mod@Mu
      sigma1 <- init.mod@Sigma
      sigma0 <- init.mod@sigma0
      e <- init.mod@e
      p <- init.mod@P[,1]
      lambda <- invlogit(1 - init.mod@zeta)
      b <- init.mod@b
      ProbMat <- init.mod@Theta
      V <- init.mod@V
      ## Sigma1 <- matrix(rep(sigma1, I), nrow = N)
      ## Sigma0 <- matrix(rep(sigma0, I), nrow = N)
      E <- matrix(rep(e, I), nrow = N)
      ## Mu1 <- matrix(rep(mu1, I), nrow = N)
      B <- matrix(rep(b, each = K), nrow = K)
      if(J > ncol(init.mod@Z)) {
        addJ <- J - ncol(init.mod@Z)
        Z <- init.mod@Z
        probz <- init.mod@probz
        while(addJ > 0) {
          addJ <- addJ - 1
          colid <- which.max(probz)
          weiMat <- matrix(runif(I * 2), nrow = I)
          weiMat <- weiMat / apply(weiMat, 1, sum)
          Z <- cbind(Z, Z[ , colid ] * weiMat[ , 2 ])
          Z[ , colid ] <- Z[ , colid ] * weiMat[ , 1 ]
          probz <- apply(Z, 2, mean)
        }
        alllik <- init.mod@alllik
        totallik <- max(init.mod@alllik)
        b.prob <- b
        predZ <- Z
        W <- init.mod@Theta[ , seq_len(I) ] %*% (Z * (1 - b)) / rep(apply(Z * (1 - b), 2, sum), each = K)

      } else {
        if(init.mod@converged & prod(struct == init.mod@Struct) == 1)
          return(init.mod) 
        W <- init.mod@W[ , 1:J ]
        predZ <- Z <- init.mod@Z[, 1:J ]
        probz <- init.mod@probz[ 1:J ]
        b.prob <- b <- init.mod@b
        alllik <- init.mod@alllik
        totallik <- max(init.mod@alllik)
      }
      clustOrder <- .orderCluster(W, struct)
      W <- W[ , clustOrder ]
      Z <- Z[ , clustOrder ]
      W <- .structure(rbind(W, 1-W), struct)[ 1:K, ]
      probz <- probz[ clustOrder ]
    }

    ## EM algorithm
    oldpar <- 1
    newpar <- c(c(W), probz, lambda, p, c(mu1), e, c(sigma1), sigma0)
    
    write.out(out, "finish initializing clusters. ")
    
    while(outitr < maxitr) {

      outitr <- outitr + 1
      if(max(abs(oldpar - newpar)) < tol) {
        break
      }

      ## transform everything into matrices
      Sigma0 <- matrix(rep(sigma0, I), nrow = N)
      E <- matrix(rep(e, I), nrow = N)
      B <- matrix(rep(b, each = K), nrow = K)

      PDF <- matrix(0, nrow = N, ncol = I * (nsig + 1))
      PDF[ , seq_len(I) ] <- logdensity(Y, Mu0 * E, Sigma0, family = family)
      for(s in seq_len(nsig)) {
        idx <- s * I + seq_len(I)
        Sigma1 <- matrix(rep(sigma1[ , s ], I), nrow = N)
        Mu1 <- matrix(rep(mu1[ , s ], I), nrow = N)
        if(family == "negbin") {
          PDF[ , idx ] <- logdensity(Y - min.count , Mu1 - min.count , Sigma1, family = family)
        } else {
          PDF[ , idx ] <- logdensity(Y, Mu1, Sigma1, family = family)
        }
      }
      PDF <- trimLogValue(PDF)
      PDF <- crossprod(t(D), PDF)

      oldlik <- totallik
      totallik <- .Call("loglik_mix", W, p, invlogit(-lambda), probz, PDF, V, package="MBASIC")
      write.out(out, paste("itr", outitr, "lik", round(totallik, 2)))

      alllik <- c(alllik, totallik)
      
      if(maxlik < totallik) {
        maxlik <- totallik
        bestb <- b.prob
        bestW <- W
        bestTheta <- ProbMat
        bestV <- V
        bestzeta <- invlogit(- lambda)
      }

      # E step
      # M step for some parameters
      PDF[ PDF < -100000 ] <- -100000
      mcmc.result <- .Call("e_step_mix", W, p, invlogit(- lambda), probz, PDF, V, package = "MBASIC")
      
      ## Expected Theta matrix
      ProbMat <- mcmc.result[["Theta_mean"]]
      ## Maximizers
      lambda <- logit(1 - mcmc.result[["zeta"]])
      p <- mcmc.result[["p"]]
      W <- mcmc.result[["W"]]
      probz <- mcmc.result[["probz"]]
      V <- mcmc.result[["V"]]
      predZ <- mcmc.result[["predZ"]]
      b.prob <- mcmc.result[["b_prob"]]

      W <- cbind(W, 1 - W)
      clustOrder <- .orderCluster(W, struct)
      oldW <- W
      for(s in 1:(ncol(W)/J))
        W[ , 1:J + (s-1) * J ] <- oldW[ , clustOrder + J * (s-1) ]
      W <- .structure(rbind(W[ , 1:J ], W[, J + (1:J) ]), struct)[1:K, ]
      probz <- probz[ clustOrder ]
      predZ <- predZ[, clustOrder ]

      ## M-step
      if(family != "lognormal") {

        for(s in seq_len(nsig)) {

          ## estimate for mu1
          ProbTheta <- crossprod(D, ProbMat[ , seq_len(I) + I * (s) ])
          F1 <- apply((Y - min.count) * ProbTheta, 1, sum)
          F2 <- apply(ProbTheta, 1, sum)
          m1 <-  F1 / F2
          m1[ m1 < 0 ] <- 0
          mu1[ , s ] <- m1 + min.count

          ## estimate sigma1

          F3 <- apply((Y - min.count) * (Y - min.count) * ProbTheta, 1, sum)
          vari <- F3 / F2 - m1 * m1
          sigma1[ ,s ] <- m1 / (vari / m1 - 1)
          sigma1[ sigma1[ ,s ] <= 0, s ] <- 100

        }

        ## estimate e
        ProbTheta <- crossprod(D, ProbMat[ , seq_len(I) ])
        F1 <-  apply(Y * ProbTheta, 1, sum)
        F2 <- apply(Mu0 * ProbTheta, 1, sum)
        e <- F1 / F2
        e[ e>1 ] <- 1
        e[ e<0.1 ] <- 0.1

        ## estimate sigma0

        E <- matrix(rep(e, I), nrow = N)
        F1 <- apply(Y * Y * ProbTheta, 1, sum)
        F2 <- apply(Mu0 * Mu0 * E * E * ProbTheta, 1, sum)
        F3 <- apply(Mu0 * E  * ProbTheta, 1, sum)
        sigma0 <- F2 / (F1 - F2 - F3)
        sigma0[ sigma0 < 0 ] <- 100

      } else {

        for(s in seq_len(nsig)) {
          ## estimate for mu1
          ProbTheta <- crossprod(D, ProbMat[ , I * (s) + seq_len(I) ])
          F1 <- apply(log(Y + 1) * ProbTheta, 1, sum)
          F2 <- apply(ProbTheta, 1, sum)
          mu1[ , s] <-  F1 / F2

          ## estimate sigma1

          sigma1[ , s] <- apply((log(Y + 1) - matrix(rep(mu1[ , s], I), nrow = N)) ^ 2 *
                                ProbTheta, 1, sum) / apply(ProbTheta, 1, sum)

          mu1[ is.na(mu1[,s ]), s ] <- mean(mu1[,s], na.rm = TRUE)
          sigma1[ is.na(sigma1[,s ]), s ] <- mean(sigma1[,s], na.rm = TRUE)

        }
        ## estimate e
        ProbTheta <- crossprod(D, ProbMat[ , seq_len(I)])
        F1 <-  apply(log(Y + 1) * ProbTheta * Mu0, 1, sum)
        F2 <- apply(Mu0 * Mu0 * ProbTheta, 1, sum)
        e <- F1 / F2
        e[ e>1 ] <- 1
        e[ e<0.1 ] <- 0.1
        e[ is.na(e) ] <- mean(e, na.rm = TRUE)

        ## estimate sigma0

        sigma0 <- apply((log(Y + 1) - Mu0 * matrix(rep(e, I), nrow = N)) ^ 2 *
                        ProbTheta, 1, sum) / apply(ProbTheta, 1, sum)
        sigma0[ is.na(sigma0) ] <- mean(sigma0, na.rm = TRUE)

      }

      lower.b <- apply(Mu0 * e, 1, mean)
      if(family == "negbin") {
          lower.b <- lower.b + 1
      }
      if(family != "negbin")
          for(s in seq_len(nsig))
              mu1[ mu1[ ,s ] < lower.b, s ] <- lower.b[ mu1[,s] < lower.b ]
      for(s in seq_len(nsig)[-1])
          mu1[ mu1[ , s ] == mu1[ , s - 1 ] , s ] <- mu1[ mu1[ , s ] == mu1[ , s - 1 ] , s - 1 ] + 0.1
   
      ## save intermediate results

      if(FALSE) {
        if((outitr %% 10 == 0) & (!is.null(outfile))) {
          tmpfit <-  list(
                          W = bestW,
                          V = bestV,
                          b = bestb,
                          Z = predZ,
                          Theta = bestTheta,
                          zeta = bestzeta,
                          aic = - 2 * tail(alllik, 1) + 2 * numpar,
                          bic = - 2 * tail(alllik, 1) + log(N * I) * numpar,
                          aicc = -2 * tail(alllik, 1) + 2 * numpar + 2 * numpar * (numpar + 1) / (N * I - numpar - 1),
                          alllik = alllik,
                          e = e,
                          mu1 = mu1,
                          sigma1 = sigma1,
                          sigma0 = sigma0,
                          p = p,
                          probz = probz
                          )
          save(tmpfit, file =outfile)
        }
      }

      oldpar <- newpar
      newpar <- c(c(W), probz, lambda, p, mu1, e, sigma1, sigma0)

    }## finish outer loop

    conv <- FALSE
    if(outitr < maxitr)
      conv <- TRUE

    rownames(bestW) <- rownames(bestTheta) <- rownames(bestV) <- unique(fac)
    
    new("MBASICFit",
        Theta = bestTheta,
        W = bestW,
        V = bestV,
        Z = predZ,
        b = bestb,
        aic = - 2 * tail(alllik, 1) + 2 * numpar,
        bic = - 2 * tail(alllik, 1) + log(N * I) * numpar,
        aicc = -2 * tail(alllik, 1) + 2 * numpar + 2 * numpar * (numpar + 1) / (N * I - numpar - 1),
        lik = tail(alllik, 1),
        alllik = alllik,
        zeta = bestzeta,
        Mu = mu1,
        Sigma = sigma1,
        sigma0 = sigma0,
        e = e,
        probz = probz,
        P = cbind(p, 1 - p),
        converged = conv,
        Struct = struct##,
##        Theta.err = NULL,
##        ARI = NULL,
##        W.err = NULL,
##        MisClassRate = NULL
      )
  }

#' @name MBASIC.binary.all
#' @title Bayesian clustering model for binary state matrix with prior estimated background means.
#'
#' @param Y An N by I matrix containing the data from N experiments across I observation units (loci).
#' @param Mu0 An N by I matrix for the prior estimated mean for the background state, for N experiments across the I observation units (loci).
#' @param fac A vector of length N denoting the experimental condition for each replicate.
#' @param allStruct A list of matrices indicating the levels of the signal matrix.
#' @param allJ A list for the numbers of clusters for each candidate model.
#' @param parent.id A vector for the identifier of the parent model for each candidate model. The parent model has a more generalized structure compared to the child model, therefore, its likelihood should be smaller.
#' @param family The distribution of family to be used. Either "lognormal" or "negbin". See details for more information.
#' @param burnin An integer value for the number of iterations in initialization. Default: 20.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for error in checking the E-M algorithm's convergence. Default: 1e-04.
#' @param nsig The number of mixture components for the distribution of the signal state.
#' @param zeta The initial value of the proportion of unclustered units. Default: 0.2.
#' @param outdir The file directory for writing the intermediate results every 10 E-M iterations. This can be useful when the running time until final convergence is long. Default: NULL (no intermediate result is saved).
#' @param ncores The number of parallel sessions used.
#' @details
#' TODO
#' @useDynLib MBASIC
#' @return A 'MBASICFit' class object.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' \dontrun{
#' }
#' @export
MBASIC.binary.fitall <- function(Y, Mu0, fac, allJ = NULL, allStruct = NULL,
                                 zeta = 0.2, maxitr = 1000, burnin = 100, family = "negbin",
                                 tol = 1E-4, nsig, min.count = 1, parent.id,
                                 outdir, ncores = 1) {
  require(ggplot2)
  require(doMC)
  registerDoMC(ncores)
  
  if(!file.exists(outdir)) {
    dir.create(outdir)
  }

  ## Assignment in Global environment because this is to be exported to clusters
  allliks <<- NULL
  if(!is.null(allJ)) {
    modelIds <- seq_along(allJ)
  } else {
    modelIds <- seq(length(allStruct))
  }

  if(length(parent.id) != length(modelIds)) {
    stop("Error: the length of 'parent.id' must be the same as the length of 'allJ' and 'allStruct'.")
  }

  if(length(parent.id) == 0) {
    stop("Error: must specify either 'allJ' or 'allStruct'.")
  }

  fitOneModel <- function(id) {
    outfile <- file.path(outdir,
                         paste("Model_", id, ".Rda", sep = ""))
    out <- file.path(outdir,
                     paste("Model_", id, ".out", sep = ""))
    thisJ <- thisStruct <- NULL
    if(!is.null(allJ)) {
      thisJ <- allJ[id]
    }
    if(!is.null(allStruct)) {
      thisStruct <- allStruct[id]
    }
    
    ## determine the parent model
    parentModel <- NULL
    parentId <- id
    if(!is.null(allliks)) {
      while(parent.id[parentId] > 0) {
        if(allliks[parent.id[parentId]] >
           allliks[parentId]) {
          parentId <- parent.id[parentId]
        } else {
          break
        }
      }
    }
    parentModelFitFile <-
      file.path(outdir,
                paste("Model_",
                      parentId,
                      ".Rda",
                      sep = ""))
    if(file.exists(parentModelFitFile)) {
      fit <- NULL
      load(parentModelFitFile)
      if(!is.null(fit)) {
        write.out(out, paste("Initialized from model", parentId, sep = ""))
        parentModel <- fit
      }
    }
    ## Fit the model
    fit <- MBASIC.binary(Y = Y,
                         Mu0 = Mu0,
                         fac = fac,
                         J = thisJ,
                         struct = thisStruct,
                         zeta = zeta,
                         maxitr = 50,
                         burnin = 20,
                         outfile = outfile,
                         out = out,
                         init.mod = parentModel,
                         family = family,
                         tol = tol,
                         nsig = nsig,
                         min.count = min.count)
    save(fit, file = outfile)
    gc()
    return(c(tail(fit@alllik, 1), fit@bic))
  }
  
  for(i in seq(as.integer(maxitr / 50))) {
    message("Pass ", i)
    if(FALSE) {
      clusterExport(cl, c("Y", "Mu0", "fac", "allJ", "allStruct", "zeta", "maxitr", "burnin", "family", "tol", "nsig", "min.count", "parent.id", "outdir", "modelIds", "allliks"))
      clusterEvalQ(cl, {library(MBASIC)})
      results <- clusterApplyLB(cl, as.list(modelIds), fitOneModel)
    }

    results <- foreach(i=modelIds) %dopar% fitOneModel(i)
    message("Pass ", i, " finished.")
    results <- matrix(unlist(results), nrow = 2)
    allliks <- results[1, ]
    allbic <- results[2, ]
    message(paste(round(allliks), collapse = " , "))
    dat <- data.frame(t(results))
    names(dat) <- c("Log_Likelihood", "BIC")
    dat$Model_ID <- seq(nrow(dat))
    ## monitor BIC and likelihood
    localenv <- environment()
    pdf(file.path(outdir, "FitScore.pdf"))
    print(ggplot(data = dat, aes(x = Model_ID, y = Log_Likelihood),
                 environment = localenv) + geom_line())
    print(ggplot(data = dat, aes(x = Model_ID, y = BIC),
                 environment = localenv) + geom_line())
    dev.off()
    gc()
  }
}

