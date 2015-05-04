#' @name MBASIC.sim.binary
#' @title Simulate data for the MBASIC model with binary states.
#' @description This function simulates MBASIC model with two states: background or binding. This mimics the standard ChIP-seq experiment.
#' @param I An integer for the total number of units.
#' @param fac A vector of length N denoting the experimental condition for each replicate.
#' @param J An integer for the number of clusters.
#' @param struct An K by J integer matrix. The j-th column denotes the levels for the cluster level parameter. See details for more information. Default: NULL.
#' @param beta Hyper parameter for simulating the entries of the W matrix. Each entry in the W matrix follows distribution Beta(beta, beta). Default: 0.1.
#' @param zeta The probability that each unit does not belong to any cluster. Default: 0.1.
#' @param beta.non The hyper parameter for simulating the probability for unclustered units. Default: 0.1.
#' @param xi Parameter for the magnitude of each observations. Default: 6. See \link{MBASIC.sim} for more information.
#' @param f Numeric value for the target fold change. Default: 5. See \link{MBASIC.sim} for more information.
#' @param family A parameter for the family of distribution to be used. Either "lognormal" or "negbin". Default: "lognormal". 
#' @return A list containing:
#' \tabular{ll}{
#'  Mu0 \tab An N by I matrix. The (n,i)-th entry is the mean of the control experiment for the n-th experiment at unit i. \cr
#'  Y \tab A matrix where each row is the counts for each replicates at all loci. \cr
#'  W \tab A K by J matrix. Each row is the indicators of the loci sets related an individual replicate. \cr
#'  Z \tab An I by J matrix. Each column is the indicator for an individual loci set.\cr
#'  Theta \tab A K by I matrix. The (k,i)-th element is the indicator of whether the i-th unit is binding for condition k.\cr
#;  e \tab A vector of length N. The n-th element is the normalization parameter for the n-th experiment.\cr
#'  pi \tab A vector of length J. The j-th entry is the prior probability for each loci to belong to the j-th cluster.\cr
#'  non.id \tab A vector of length I indicating the number of each column not to belong to any cluster.\cr
#' fac \tab A vector of length N denoting the experimental condition for each replicate.\cr
#'  bkng \tab Mean value of the responses with background state.\cr
#'  snr \tab The ratio between the mean value of the responses with the binding state and the background state.\cr
#' }
#' @seealso \code{\link{MBASIC.sim}}
#' @examples
#' dat.sim <- MBASIC.sim.binary(I = 100, fac = rep(1:5, each = 2), J = 3, f = 5)
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @importFrom gtools rdirichlet
#' @importFrom msm rtnorm
#' @export
MBASIC.sim.binary <- function(I, fac, J, f, struct=NULL, beta=0.1, zeta = 0.1, beta.non=0.1, xi = 6, family="lognormal") {
  
        ## prespecified
        K <- length(unique(fac))
        if(!is.null(struct)) {
          J <- ncol(struct)
          if(nrow(struct) != K) {
            stop("Error: dimension of 'struct' is not consistent with 'K'.")
          }
        } else {
          struct <- matrix(1:K, nrow = K, ncol = J)
        }

        ##prior parameters
        alpha <- 200##dirichlet prior
        alpha0 <- 400
        omega <- 10
        nu <- 10 / (omega - 1)
        tau <-  0.1
        eta <- 0.75
        epsilon <- 0.1 ^ 2
        N <- length(fac)
        
        ## design matrix Z is K by sum(n)
        D <- matrix(0, nrow = K, ncol = N)
        for(k in 1:K) {
                D[ k, fac == unique(fac)[ k ] ] <- 1
        }
        
        delta <- rdirichlet(1, rep(alpha, J))
        
        ## Initialize parameters
        Z <- matrix(0, ncol = J , nrow = I) 
        for(i in 1:I) {
                Z[ i, ] <- rmultinom(1, size = 1, prob = delta)
        }

        ## allset <- max(apply(Z, 2, sum))
        W <- matrix(0, nrow = K, ncol = J)
        for(j in 1:J) {
          if(length(unique(struct[,j])) == 2) {
            bindp <- rbeta(2, beta, beta)
            lev <- unique(struct[,j])
            if(j %% 2 == 1) {
              W[ struct[,j] == lev[1], j ] <- max(c(bindp[1], 1-bindp[1]))
              W[ struct[,j] == lev[2], j ] <- min(c(bindp[2], 1-bindp[2]))
            }
            else {
              W[ struct[,j] == lev[1], j ] <- min(c(bindp[1], 1-bindp[1]))
              W[ struct[,j] == lev[2], j ] <- max(c(bindp[2], 1-bindp[2]))
            }
          } else{
            for(l in unique(struct[,j])) {
              W[ struct[,j] == l, j ] <- rbeta(1, beta, beta)
            }
          }
        }

        ProbTheta <- tcrossprod(W, Z)

        b <- rbinom(I, size = 1, prob = zeta)
        non.id <- which(b == 1)
        
        ProbTheta[ , b == 1 ] <- rep(rbeta(length(non.id), beta.non, beta.non), each = K)

        Theta <- matrix(rbinom(K * I, size = 1, prob= ProbTheta), nrow = K, ncol = I)
        
        sigma1 <- 1 / rgamma(N, shape = omega, scale = nu)
        sigma0 <- 1 / rgamma(N, shape = omega, scale = nu)
        size0 <- runif(N, 5, 10)
        size1 <- runif(N, 5, 10)
        
        e <- rtnorm(N, mean = eta, sd = sqrt(epsilon), lower = 0.5, upper = 1)

        Mu <- Mu0 <- matrix(0, nrow = N, ncol = I)
        Delta <- crossprod(D, Theta)
        
        Sigma1 <- matrix(rep(sigma1, I), nrow = N)
        Sigma0 <- matrix(rep(sigma0, I), nrow = N)
        Sigma <- Sigma0
        Size <- matrix(0, nrow = N, ncol = I)

        .var2size <- function(mu, lnsd) {
          mu <- exp(mu)
          v <- mu * mu * (exp(lnsd * lnsd * 2) - exp(lnsd * lnsd))
          if(v < mu)
            return(1000)
          return(mu / (v / mu - 1))
        }
        
        for(l in 1:N) {
          Mu[ l, ] <- Mu0[ l, ] <- rtnorm(I, mean = xi, sd = sqrt(tau), lower = 0, upper = xi * 2)
          if(family == "lognormal") {
            Mu[ l, ] <- Mu0[ l, ] * e[ l ]
            Mu[ l, Delta[ l, ] == 1 ] <- (mean(Mu[l,]) + log(f))
          } else{
            Mu[ l, ] <- Mu0[ l, ] * e[ l ]
            Mu[ l, Delta[ l, ] == 1 ] <- (mean(Mu[l,]) * f)
          }

          Size[ l, Delta[ l, ] == 1 ] <- .var2size(mean(Mu[ l, Delta[ l, ] == 1 ]), sigma1[ l ])
          Size[ l, Delta[ l, ] == 0 ] <- .var2size(mean(Mu[ l, Delta[ l, ] == 0 ]), sigma0[ l ])
        }

        if(family == "lognormal")
          Y <-  matrix(
                                        #          rnbinom(N * I, mu = as.vector(Mu), size = 5),
                       as.integer(
                                  exp(
                                      rtnorm(N * I, mean = as.vector(Mu), sd = Sigma, lower = 0)
                                     )
                                 ) - 1,
                       nrow = N
                      )
        else
          Y <- matrix(
                      rnbinom(N * I, mu = Mu, size = Size),
                      nrow = N)


        if(FALSE) {
          E <- matrix(rep(e, I), nrow = N)
          Mu1 <- matrix(rep(mu1, I), nrow = N)
          B <- matrix(rep(b, each = K), nrow = K)
          pW <- tcrossprod(W, delta) * (1 - zeta) + zeta * p.non  
          pW <- rep(pW, n)
          totallik <- sum(log(rep(pW, I) * exp(logdnorm(log(Y + 1), Mu1, Sigma1)) + rep(1-pW, I) * exp(logdnorm(log(Y + 1), Mu0 * E, Sigma0))))
        }
       
	bkng <- mean(Y[ Delta == 0 ])
	snr <- mean(Y[ Delta == 1 ]) / bkng

        return(list(Mu0 = Mu0, Y = Y, W = W, Z = Z, pi = delta, Theta = Theta, e = e, non.id=non.id, fac = fac, bkng = bkng, snr = snr))
        
}

#' Simulate for the Theta matrix
#' @name MBASIC.sim.state
#' @title Simulate the state matrix for general MBASIC models.
#' @param I An integer for the total number of units.
#' @param K An integer for the number of different experimental conditions.
#' @param J An integer for the number of clusters.
#' @param S An integer for the number of states.
#' @param struct An K by J integer matrix. The j-th column denotes the levels for the cluster level parameter. See details for more information. Default: NULL.
#' @param statemap A vector consisted of 1, 2, ..., S, representing which state each component belongs to. Default: NULL.
#' @param delta A vector the same length as statemap, or NULL. This is the dirichlet prior parameter to simulate the probability across the components for each CLUSTERED unit and each experiment. If NULL, rep(0.1, length(statemap)) is used.
#' @param delta.non  A vector of length S, or NULL. This is the dirichlet prior parameter to simulate the probability across the components for each UNCLUSTERED unit and each experiment. If NULL, rep(0.1, length(statemap)) is used.
#' @param zeta The probability that each unit does not belong to any cluster. Default: 0.1.
#' @return A list containing:
#' \tabular{ll}{
#'  Theta \tab A K by I matrix. The (k,i)-th element is the indicator of whether the i-th unit is binding for condition k.\cr
#'  W \tab A K by J matrix. Each row is the indicators of the loci sets related an individual replicate. \cr
#'  Z \tab An I by J matrix. Each column is the indicator for an individual loci set.\cr
#'  delta \tab Same as the input argument \code{delta}.\cr
#'  zeta \tab Same as the input argumnent \code{zeta}.\cr
#'  non.id \tab A vector of length I indicating the number of each unit not to belong to any cluster.\cr
#' }
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' dat.sim <- MBASIC.sim.state(I = 100, K = 10, J = 3)
#' @importFrom gtools rdirichlet
#' @export
MBASIC.sim.state <- function(I, K, J, S = 2, struct = NULL, delta = NULL, delta.non = NULL, zeta = 0.1, statemap = NULL) {
  
  if(is.null(statemap)) {
    statemap <- seq(S)
  } else if(prod(sort(unique(statemap)) == seq(S)) != 1) {
    stop("Error: statemap must be a vector consisted of 1, 2, ..., S.")
  }

  ##prior parameters
  if(is.null(delta))
    delta <- rep(0.1, length(statemap))
  else if(length(statemap) != length(delta)) {
  stop("delta ", paste(delta, collapse = " ")," statemap ", paste(statemap, collapse = " "), " Error: length of 'delta' must be the same as 'statemap'.")
  }

  if(is.null(delta.non))
    delta.non <- rep(0.1, length(statemap))
  else if(length(statemap) != length(delta.non)) {
  stop("delta.non ", paste(delta.non, collapse = " ")," statemap ", paste(statemap, collapse = " "), " error: length of 'delta.non' must be the same as 'statemap'.")
  }

  if(!is.null(struct))
    if(ncol(struct) != J | nrow(struct) != K)
      stop("Error: matrix 'struct' is not of the correct dimension.")

  stateMap <- matrix(0, nrow = length(statemap), ncol = S)
  stateMap[cbind(seq_along(statemap), statemap)] <- 1

  delta <- c(t(stateMap) %*% delta)
  delta.non <- c(t(stateMap) %*% delta.non)
  
  W <- matrix(t(rdirichlet(K * J, delta)), nrow = S * K, ncol = J)

  for(j in 1:J)
    if(length(unique(struct[, j ])) < K)
      for(l in unique(struct[ , j ])) {
        idx <- rep((which(struct[ , j ] == l) - 1) * S, each = S) + seq_len(S)
        W[ idx, j ] <- W[ idx[ seq_len(S) ], j ]
      }
  
  ## Initialize parameters
  Z <- matrix(0, ncol = J, nrow = I) 
  for(i in 1:I) {
    Z[ i, ] <- rmultinom(1, size = 1, prob = rep(1, J))
  }

  ProbMat <- tcrossprod(W, Z)

  non.id <- which(rbinom(I, size=1, prob = zeta) == 1)
  
  ProbMat[ , non.id ] <- c(t(rdirichlet(K * length(non.id), delta.non)))
  
  Theta <- matrix(apply(matrix(ProbMat, nrow = S), 2, function(x) which(rmultinom(1, size = 1, prob = x) == 1)), nrow = K)

  return(list(Theta = Theta, W = W, Z = Z, delta = delta, zeta = zeta, non.id = non.id))
        
}

#' @name MBASIC.sim
#' @title Simulate data for the general MBASIC model.
#' @param xi Parameter for the magnitude of each observations. See details for more information.
#' @param family A parameter for the family of distribution to be used, must be "lognormal", "negbin" or "binom". Default: "lognormal".
#' @param f A numerical value that determine the difference of the means between different states. See details for more information. Default: 5.
#' @param I An integer for the total number of units.
#' @param fac A vector of length N denoting the experimental condition for each replicate.
#' @param J An integer for the number of clusters.
#' @param S An integer for the number of states. Default: 2.
#' @param struct An K by J integer matrix. The j-th column denotes the levels for the cluster level parameter. See details for more information. Default: NULL.
#' @param delta A vector of length S, or NULL. This is the dirichlet prior parameter to simulate the probability across the S states for each CLUSTERED unit and each experiment. If NULL, rep(0.1,S) is used.
#' @param delta.non  A vector of length S, or NULL. This is the dirichlet prior parameter to simulate the probability across the S states for each UNCLUSTERED unit and each experiment. If NULL, rep(0.1,S) is used.
#' @param zeta The probability that each unit does not belong to any cluster. Default: 0.1.
#' @details
#' MBASIC.sim allows three types of distributions:\cr
#' For the "lognormal" family, entries in the matrix Y follows distribution: log(Y[n,i] + 1) | Theta[n,i]=s ~ N(Mu[n,s], stdev[s]).\cr
#' For the "negbin" family, entries in the matrix Y follows distribution: Y[n,i] | Theta[n,i]=s ~ NB(Mu[n,s], stdev[s]).\cr
#' For the "binom" family, entries in the matrices X and Y follows distribution: X[n,i] ~ Pois(xi), Y[n,i]| Theta[n,i]=s,X[n,i] ~ Binom(X[n,i],Mu[n,s]). In this package, NB(mu,size) denotes a Negative binomial distribution with mean mu and variance mu(1+mu/size).\cr
#' For "lognormal" or "negbin" families, Mu[n,s]~N(prior.mean[s],prior.sd[s]). Hyper paramters prior.mean and prior.sd are set differently under the two distributional families. For the "lognormal" family, where prior.mean[s] = xi+log((s-1)(f-1)+1), and prior.sd=log(f)/30. For the "negbin" family, prior.mean[s]=xi*((s-1)(f-1)+1), and prior.sd=(f-1)*xi/6. In general, xi is the mean for the state S=1, and f is roughly the ratio between the means from state S=2 and S=1.\cr
#' For the "binom" family, Mu[n,s] ～ Beta(s * f, (S + 1 - s) * f).
#' @return A list containing:
#' \tabular{ll}{
#'  Y \tab An N by I matrix. The (n,i)-th entry is the observed value at the i-th unit for the n-th experiment. \cr
#'  X \tab An N by I matrix for the "binom" family, or NULL otherwise. The (n,i)-th entry is the observed size parameter for the i-th unit in the n-th experiment. \cr
#'  fac \tab Same as the input argument \code{fac}.\cr
#'  Theta \tab A K by I matrix. The (k,i)-th element is the indicator of whether the i-th unit is binding for condition k.\cr
#'  W \tab A K by (J*S) matrix. The (k,J*(s-1)+j)-th entry is the probability of units in the j-th cluster have state s under the k-th experimental condition. \cr
#'  Z \tab An I by J matrix. Each column is the indicator for an individual loci set.\cr
#'  delta \tab Same as the input argument \code{delta}.\cr
#'  zeta \tab Same as the input argumnent \code{zeta}.\cr
#'  non.id \tab A vector of length I indicating the number of each unit not to belong to any cluster.\cr
#'  prior.mean \tab A vector of length S, the hyper means for the S states for each experiment.\cr
#'  prior.sd \tab A vector of length S, the hyper parameters for the dispersions for the S states for each experiment.\cr
#'  Mu \tab A N by S matrix. The (n,s)-th entry is the mean values of the response for the s-th state in the n-th experiment.\cr
#'  stdev \tab A vector of length S. The dispersion parameter for the S states for the observed data.\cr
#' }
#' @seealso \code{\link{MBASIC.sim.state}}
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' dat.sim <- MBASIC.sim(xi = 2, I = 100, fac = rep(1:5, each = 2), J = 3)
#' @importFrom gtools rdirichlet
#' @importFrom msm rtnorm
#' @export
MBASIC.sim<- function(xi, family = "lognormal", struct = NULL, I, fac, J, S = 2, f = 5, delta = NULL, delta.non = NULL, zeta = 0.1, statemap = NULL) {
  
  if(!family %in% c("lognormal", "negbin", "binom")) {
    stop("Error: 'family' must be one of 'lognormal', 'negbin' or 'binom'.")
  }

  if(is.null(statemap)) {
    statemap <- seq(S)
  }
 
  if(prod(sort(unique(statemap)) == seq(S)) != 1) {
    stop("Error: statemap must be consisted of 1, 2, ... S.")
  }

  M <- length(statemap)
  K <- length(unique(fac))
  N <- length(fac)

  para.theta <- MBASIC.sim.state(I=I, K=K, J=J, S=S, delta=delta, delta.non=delta.non, zeta=zeta, struct = struct, statemap = statemap)
   V <- matrix(0, nrow = N, ncol = M)
  for(s in seq(S)) {
    ids <- which(statemap == s)
    if(length(ids) == 1)
      V[, ids] <- 1
    else
      V[, ids] <- t(matrix(rdirichlet(N, delta[ids]), nrow = length(ids)))
  }
  
  Delta <- matrix(0, nrow = N, ncol = I)
  for(n in seq(N)) {
    for(s in seq(S)) {
      ids <- which(para.theta$Theta[fac[n], ] == s)
      if(length(ids) == 0)
        next
      if(sum(statemap == s) == 1) {
        Delta[n, ids] <- which(statemap == s)
      } else {
        Delta[n, ids] <- sample(which(statemap == s), length(ids), replace = TRUE, prob = V[n, statemap == s])
      }
    }
  }
  
  if(family == "lognormal") {
    prior_mean <- xi + log ((seq(M) - 1) * (f - 1) + 1)
    prior_sd <- 0.05
    ## sdev <- diff(exp(prior_mean))[1]
    ## stdev <- sqrt((log(sdev ^ 2 + exp(2 * prior_mean)) - 2 * prior_mean) / 2)
    stdev <- 0.5
  } else {
    prior_mean <- xi * ((seq(M) - 1) * (f - 1) + 1)
    prior_sd <- 0.5
    sdev <- prior_mean + sqrt(prior_mean)
    stdev <- prior_mean / (sdev / prior_mean - 1)
    stdev[stdev < 0] <- 5
    stdev[stdev > 5] <- 5
  }
    
  if(family != "binom") {  
    prior.Mu <- t(matrix(rtnorm(M * N, mean = prior_mean, sd = prior_sd, lower = 0), nrow = M))
  } else {
    prior.Mu <- t(matrix(rbeta(M * N, f * seq(M), f * rev(seq(M))), nrow = M))
  }

  ## order the means within each replicate
  for(n in seq(N)) {
    od <- order(prior.Mu[n, ])
    prior.Mu[n, ] <- prior.Mu[n, od]
  }
  
  Mu <- matrix(0, nrow = N, ncol = I)
  for(n in 1:N)
    for(m in seq(M))
      Mu[ n, Delta[n,] == m ] <- prior.Mu[ n, m ]

  X <- NULL
  if(family == "lognormal") {
    Y <- matrix(as.integer(exp(rtnorm(N * I, mean = Mu, sd = stdev, upper = max(Mu) + 1))), nrow = N)
  } else if(family == "negbin") {
    Y <- matrix(rnbinom(N * I, mu = Mu, size = stdev), nrow = N)
  } else {
    X <- matrix(rpois(N * I, lambda = xi), nrow = N)
    Y <- matrix(rbinom(N * I, size = X, prob = Mu), nrow = N)
  }
  
  bkng <- mean(Y[ Delta == 1 ])
  snr <- mean(Y[ Delta == 2 ]) / mean(Y[ Delta == 1 ])

  return(list(Theta = para.theta$Theta, Y = Y, X = X, fac = fac, W = para.theta$W, Z = para.theta$Z, V = V, delta = para.theta$delta, zeta = para.theta$zeta, prior.mean = prior_mean, prior.sd = prior_sd, stdev = stdev, Mu = prior.Mu, bkng = bkng, snr = snr, non.id = para.theta$non.id))
  
}
