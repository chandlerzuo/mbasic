#' @name MBASIC
#' @title Bayesian clustering model for a state-space matrix.
#'
#' @description This function is designed to analyze general state-space models. The data consists of observations over I units under N experiments with K different conditions. There are S states for each experiment and unit.
#' @param Y An N by I matrix containing the data from N experiments across I observation units.
#' @param S An integer for the number of states.
#' @param fac A vector of levels representing the conditions of each replicate.
#' @param struct A K by J matrix indicating the structures of each cluster.
#' @param J The number of clusters to be identified.
#' @param family The distribution of family to be used. Either "lognormal" or "negbin". See details for more information.
#' @param method A string for the fitting method, 'naive', 'em' (default) or '2em'. See details for more information.
#' @param para A list object that contains the true model parameters. Default: NULL. See details for more information.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for error in checking the E-M algorithm's convergence. Default: 1e-04.
#' @param zeta The initial value for the proportion of units that are not clustered. Default: 0.1.
#' @param out The file directory for writing fitting information in each E-M iteration. Default: NULL ( no information is outputted ).
#' @details
#' Function MBASIC currently supports two different distributional families: log-normal and negative binomial. This should be specified by the 'family' argument.\cr
#' For the log-normal distributions, log(Y+1) is modeled as normal distributions. For experiment n, if locus i has state s, distribution for log(Y[n,i]+1) is N( Mu[n,s], Sigma[n,s] ).\cr
#' For the negative binomial distributions, the meanings of Mu and Sigma are different. For experiment n, if locus i has state s, distribution of Y[n,i] is NB( Mu[n,s], Sigma[n,s] ). In this package, NB( mu, a ) denotes the negative-binomial distribution with mean mu and size a (i.e. the variance is mu*(1+mu/a) ).\cr
#'  The 'method' argument determines what fitting method will be used. The default is 'em', where the states and the clustering are simultaneously estimated. 'naive' and '2em' methods use 2-step algorithms. In the first step, both estimate the states for each unit by an E-M algorithm for each experiment. In the second step, 'naive' uses hierarchical clustering to cluster the units, while '2em' uses function 'MBASIC.state' to identify clusters.\cr
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
#' dat.sim <- MBASIC.sim( xi = 2, family = "lognormal", I = 1000, fac = rep( 1:10, each = 2 ), J = 3, S = 3, zeta = 0.1 )
#' ## Fit the model
#' dat.sim.fit <- MBASIC( Y = dat.sim$Y, S = 3, fac = rep( 1:10, each = 2), J = 3, maxitr = 3, para = NULL, family = "lognormal", method = "em", zeta = 0.1, tol = 1e-04)
#' @useDynLib MBASIC
#' @export
MBASIC <- function( Y, S, fac, J=NULL, maxitr = 100, struct = NULL, para = NULL,  family="lognormal", method = "em", zeta = 0.1, tol = 1e-4, out = NULL ){

  write.out( out, "Started" )
  
    if( ! method %in% c( "naive", "2em", "em" ) ){
        message( "Error: method must be of naive, 2em, em." )
        return
    }
    
    ## prespecified
    K <- length( unique( fac ) )
    I <- ncol( Y )
    N <- nrow( Y )
    if( length( fac ) != N )
      message( "Error: total number of replicates do not match with the number of rows in Y" )
  
  ## index
  ##    id <- rep( n[1], K )
  ##    for( k in 2:K ){
  ##        id[k] <- id[k-1] + n[k]
  ##    }
    
    ## design matrix D is K by N
    D <- matrix( 0, nrow = K, ncol = N )
    for( k in 1:K ){
      D[ k, fac == unique( fac )[ k ] ] <- 1
    }
    SampleToExp <- apply( D, 2, function( x ) which( x==1 ) )

    if( is.null( struct ) ){
        if( is.null( J ) )
            message( "Error: either struct or J must not be missing." )
        struct <- matrix( seq_len( K ), nrow = K, ncol = J )
    } else {
        if( is.null( J ) )
            J <- ncol( struct )
        J <- sum( J )
        if( ncol( struct )!= sum( J ) | nrow( struct ) != K )
            message( "Error: the dimension of struct is inconsistent with grouping structure!" )
    }
    
    numpar <- 2 * S * N + I * ( S - 1 ) + J + sum( apply( struct, 2, function(x) length( unique( x ) ) ) ) * ( S - 1 )
    
    outitr <- 0
    totallik <- oldlik <- 0
    alllik <- allerr <- allzeta <- bestW <- allmisclass <- matchId1 <- W.err <- matchId2 <- allari <- numeric( 0 )
    maxlik <- -Inf

  write.out( out, "Initialized parameters" )
  
    ## initialize mu and variance
    Sigma <- Mu <- matrix( 0, nrow = N, ncol = S )
    
    if( family == "lognormal" ){
        for( s in 1:S ){
            Y.sec <- c( Y )[ c( Y ) < quantile( c(Y), s / S ) & c( Y ) > quantile( c( Y ), ( s - 1 ) / S ) ]
            Mu[ , s] <- mean( log( Y.sec + 1 ) )
            Sigma[ , s ] <- sd( log( Y.sec + 1 ) )
        }
    } else {
        for( s in 1:S ){
            Y.sec <- c( Y )[ c( Y ) < quantile( c(Y), s / S ) & c( Y ) > quantile( c( Y ), ( s - 1 ) / S ) ]
            Mu[ , s] <- m1 <- mean( Y.sec )
            vari <- var( Y.sec )
            size <- m1 / ( vari / m1 - 1 )
            if( size > 0 )
                Sigma[ , s] <- size
            else
                Sigma[ , s ] <- 100
        }
    }

    b <- rep( 0, I )
    B <- matrix( rep( b, each = K ), nrow = K )
    
    ## initialize the matrices by hierarchical clustering
    ## in constructing Z, cluster all locis
    ## This gives deterministic initialization
    totalF <- matrix( 0, nrow = K, ncol = I )
    if( family == "lognormal" ){
        F1 <- matrix( 0, nrow = K * S, ncol = I )
        for( s in 1:S ){
            idx <- (s-1) * K + seq_len( K )
            F1[ idx, ] <-  exp( crossprod( t( D ) , - ( log( Y + 1 ) - Mu[ ,s] ) ^ 2 / 2 / Sigma[ ,s] ) )
            totalF <- totalF + F1[ idx, ]
        }
   } else {
        F1 <- matrix( 0, nrow = K * S, ncol = I )
        for( s in 1:S ){
            idx <- (s-1) * K + seq_len( K )
            F1[ idx, ] <-  exp( crossprod( t( D ), log( dnbinom( Y, size = Sigma[,s], mu = Mu[,s] ) ) ) )
            totalF <- totalF + F1[ idx, ]
        }
    }
    totalF <- t( matrix( rep( c( t( totalF ) ), S ), nrow = I ) )
    ProbMat <- F1 / totalF
    ProbMat[ totalF == 0 ] <- 1/S
    maxProb <- max( na.omit( ProbMat[ ProbMat != 1 ] ) )
    minProb <- min( na.omit( ProbMat[ ProbMat != 0 ] ) )
    ProbMat[ ProbMat > maxProb ] <- max( c( 0.999, na.omit( maxProb ) ) )
    ProbMat[ ProbMat < minProb ] <- min( c( 0.001, na.omit( minProb ) ) )
  ProbMat[ is.na( ProbMat )] <- mean( ProbMat, na.rm = TRUE )

    if( method != "em" ){
        ## naive or 2-step method
        Pi <- matrix( 1/S, nrow = K, ncol = S )
        ## em step to estimate Theta
        allpar <- c( c( Mu), c( Sigma ), c( Pi ) )
        oldpar <- 0
        for( itr in 1:maxitr ){
            ## check for convergence
            if( max( abs( oldpar - allpar ) )< tol )
                break

            ## M step
            for( s in seq_len( S ) ){
                idx <- seq_len( K ) + (s-1 ) * K
                Pi[ ,s ] <- apply( ProbMat[ idx, ], 1 ,mean )
            }

            if( family == "lognormal" ){
                for( s in seq_len( S ) ){
                    idx <- seq_len( N ) + (s-1) * N
                    Mu[ ,s ] <- apply( log( Y + 1 ) * ProbMat[ SampleToExp + ( s-1 ) * K, ], 1, sum ) / apply( ProbMat[ SampleToExp + ( s-1 ) * K, ], 1, sum )
                    M2 <- apply( log( Y + 1 ) ^ 2 * ProbMat[ SampleToExp + ( s-1 ) * K, ], 1, sum ) / apply( ProbMat[ SampleToExp + ( s-1 ) * K, ], 1, sum )
                    Sigma[ ,s ] <- M2 - Mu[ ,s ] ^2
                    Sigma[ Sigma < 0.01 ] <- 0.01
                }
            } else {
                ## negative binomial family
                for( s in seq_len( S ) ){
                    idx <- seq_len( N ) + (s-1) * N
                    Mu[ ,s ] <- apply( Y * ProbMat[ SampleToExp + ( s-1 ) * K, ], 1, sum ) / apply( ProbMat[ SampleToExp + ( s-1 ) * K, ], 1, sum )
                    M2 <- apply( Y ^ 2 * ProbMat[ SampleToExp + ( s-1 ) * K, ], 1, sum ) / apply( ProbMat[ SampleToExp + ( s-1 ) * K, ], 1, sum )
                    M2 <- M2 - Mu[ , s] ^ 2
                    Sigma[ ,s ] <- Mu[ ,s ] / ( M2 / Mu[ , s] - 1 )
                    Sigma[ Sigma[ ,s ] < 0,] <- 100
                }
            }
            ## order the means
            od <-  apply( Mu, 1, order )
            Mu <- matrix( Mu[ cbind( rep( seq_len( N ), each = S ), c( od ) ) ], ncol = S, byrow = TRUE )
            Sigma <- matrix( Sigma[ cbind( rep( seq_len( N ), each = S ), c( od ) ) ], ncol = S, byrow = TRUE )

            ## E step
            if( family == "lognormal" ){
                F1  <- matrix( 0, nrow = K * S, ncol = I )
                for( s in 1:S ){
                    idx <- (s-1) * K + seq_len( K )
                    F1[ idx, ] <-  crossprod( t( D ), - ( log( Y + 1 ) - Mu[ ,s] ) ^ 2 / 2 / Sigma[ ,s] - log( Sigma[ , s ] ) / 2 ) + log( Pi[ ,s ] )
                }
            } else {
                F1 <- matrix( 0, nrow = K * S, ncol = I )
                for( s in 1:S ){
                    idx <- (s-1) * K + seq_len( K )
                    F1[ idx, ] <-  crossprod( t( D ), dnbinom( Y, size = Sigma[,s], mu = Mu[,s], log = TRUE ) ) + log( Pi[ ,s ] )
                }
            }
            F1[ is.na( F1 ) ] <- -5000
            F1[ F1 <  -5000 ] <- -5000
            F.max <- t( matrix( apply( matrix( t( F1 ), ncol = S ), 1, max ), nrow = I ) )
            totalF <- matrix( 0, nrow = K, ncol = I )
            for( s in seq_len( S ) ){
                idx <- seq_len( K ) + ( s-1 ) * K
                F1[ idx, ] <- exp( F1[ idx, ] - F.max )
                totalF <- totalF + F1[ idx, ]
            }
            totalF <- t( matrix( rep( c( t( totalF ) ), S ), nrow = I ) )

            ProbMat <- F1 / totalF
            maxProb <- max( ProbMat[ ProbMat != 1 ] )
            minProb <- min( ProbMat[ ProbMat != 0 ] )
            ProbMat[ ProbMat > maxProb ] <- maxProb
            ProbMat[ ProbMat < minProb ] <- minProb
        }## finish iteration

        Theta <- matrix( -1, nrow = K, ncol = I )
        for( i in seq_len( K ) ){
            idx <- i + K * ( seq_len( S ) - 1 )
            Theta[ i, ] <- apply( ProbMat[ idx, ], 2, which.max )
        }
        if( !is.null( para ) )
            allerr <- mean( Theta != para$Theta )
        
        ret <- MBASIC.state( Theta, J=J, zeta = zeta, struct = struct, method = method, maxitr = maxitr, tol = tol, para = para, out = out )

        conv <- FALSE
        if( ret@converged & itr < maxitr )
          conv <- TRUE
        
        ## Pi is the proportion for components in the k experiment to have state s
        ## Pi is different from Z. Z is the posterior probability.
        
        return( new( "MBASICFit",
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
                    Theta.err = tail( allerr, 1 ),
                    ARI = ret@ARI,
                    W.err = ret@W.err,
                    MisClassRate = ret@MisClassRate
                    )
               )

    }
    
    ## initialize W, Z, b
    ## ProbMat <- D.rep %*% ProbMat
    d <- dist( t( ProbMat ) )
    mind <- apply( as.matrix(d), 1, function( x ) min( x[x>0] ) )
    thr <- quantile( mind, 1 - zeta )
    id <- which( mind < thr )
    b <- rep(1, I )
    b[id] <- 0
    d <- dist( t( ProbMat[,id] ) )
    fit <- hclust( d )
    groups <- cutree( fit, k = J )
    Z <- matrix( 0, nrow = I, ncol = J )
    Z[ cbind( 1:I, sample( 1:J, I, replace = TRUE ) )] <- 1
    Z[ id, ] <- 0
    Z[ cbind( id, groups ) ] <- 1
    W <- matrix( 1/S, nrow = K * S, ncol = J )
    for( j in seq_len( J ) ){
        if( length( id[ groups == j ] ) > 0 ){
            W[ , j ] <- apply( t( ProbMat[  , id[ groups == j ] ] ), 2, mean )
        }
    }
    predZ <- Z
    b.prob <- b
    clustOrder <- .orderCluster( W, struct )
    W <- W[ , clustOrder ]
    Z <- Z[ , clustOrder ]
    W <- .structure( W, struct )
    ## initialize p, probz
    P <- matrix( 0, nrow = I, ncol = S )
    for( s in seq_len( S ) ){
        idx <- seq_len( K ) + K * ( s - 1 )
        P[ ,s ] <- apply( ProbMat[ idx, ], 2, mean )
    }
    probz <- apply( rbind( Z, diag( rep( 1, J ) ) ), 2, mean )

    ## EM algorithm
    oldpar <- 0
    newpar <- c( c( W ), probz, zeta, c( P ), c( Mu ), c( Sigma ) )
   
    for( outitr in seq_len( maxitr ) ){
        
        if( max( abs( newpar - oldpar ) ) < tol ){
            break
        }
        
        ## transform everything into matrices
        B <- matrix( rep( b, each = K ), nrow = K )

        PDF <- matrix( 0, nrow = N, ncol = I * S )
        W.lik <- matrix( 0, nrow = K, ncol = J * S )
        for( s in seq_len( S ) ){
            PDF[ , seq_len( I ) + I * ( s - 1 ) ] <- logdensity( Y, Mu[ , s ], Sigma[ ,s ], family )
            W.lik[ , seq_len( J ) + J * ( s - 1 ) ] <- W[ seq_len( K ) + K * ( s - 1 ) , ]
        }
        PDF[ PDF > 5 ] <- 5
        PDF[ PDF < -5000 ] <- -5000
        PDF[ is.na( PDF ) ] <- mean( PDF, na.rm = TRUE )
        PDF <- crossprod( t( D ), PDF )
        
        oldlik <- totallik
        totallik <- .Call( "loglik", W.lik, P, zeta, probz, PDF, package="MBASIC" )
        write.out( out, paste( "itr", outitr, "lik", round( totallik, 2 ), "zeta", round( zeta, 2 ) ) )
        ##write.out( out, paste( "loglik in C =", totallik ) )
        
        alllik <- c( alllik, totallik )
        allzeta <- c( allzeta, zeta )
        Theta <- matrix( -1, nrow = K, ncol = I )
        for( i in seq_len( K ) ){
            idx <- i + K * ( seq_len( S ) - 1 )
            Theta[ i, ] <- apply( ProbMat[ idx, ], 2, which.max )
        }
        if( length(names(para)) != 0 ){
            allerr <- c( allerr, mean( para$Theta != Theta ) )

            ## compute misclassification rate
            W.f <- matrix( 0, nrow = K * S, ncol = J )
            for( s in seq_len( S ) )
              W.f[ s + S * seq( 0, K - 1 ), ] <- W[ seq_len( K ) + K * ( s - 1 ), ]
            
            mc <- matchCluster( W.f, para$W, predZ, para$Z, b.prob, para$non.id )
            
            write.out( out, paste( "mis-class rate ", mc$mcr ) )
            write.out( out, paste( "Error for W ",  round( mc$W.err, 3 ) ) )
            allmisclass <- c( allmisclass, mc$mcr )
            W.err <- c( W.err, mc$W.err )
            allari <- c( allari, mc$ari )
            write.out( out, paste( "ARI ", mc$ari ) )
            write.out( out, paste( "loglik", totallik, "err", round( allerr[ length( allerr ) ], 2 ) ) )

            ##matchId1 <- mc$matchId1
            ##matchId2 <- mc$matchId2
            
        }
        
        if( maxlik < totallik ){
            maxlik <- totallik
            bestb <- b.prob
            bestW <- W
            bestTheta <- Theta
        }

        ## E step
        ## M step for some parameters
        mcmc.result <- .Call( "e_step", W.lik, P, zeta, probz, PDF, package = "MBASIC" )
        
        ## Expected Theta matrix
        ProbMat <- mcmc.result[["Theta_mean"]]
        ## Maximizers
        zeta <- mcmc.result[["zeta"]]
        P <- mcmc.result[["P"]]
        W <- mcmc.result[["W"]]
        probz <- mcmc.result[["probz"]]
    
        predZ <- mcmc.result[["predZ"]]
        b.prob <- mcmc.result[["b_prob"]]

        W.aug <- matrix( 0, nrow = K * S, ncol = J )
        for( s in 1:S ){
            W.aug[ seq_len( K ) + K * ( s - 1 ), ] <- W[ ,  seq_len( J ) + J * ( s - 1 ) ]
        }
        clustOrder <- .orderCluster( W.aug, struct )
        W <- W.aug[ , clustOrder ]
        W <- .structure( W, struct )
        probz <- probz[ clustOrder ]
        predZ <- predZ[, clustOrder]
    
        ## M-step
        if( family != "lognormal" ){

            for( s in seq_len( S ) ){
                ## estimate for mu1
                idx <- seq_len( I ) + I * ( s - 1 )
                F1 <- apply( Y * crossprod( D, ProbMat[ , idx ] ), 1, sum )
                F2 <- apply( crossprod( D, ProbMat[ , idx ] ), 1, sum )
                ## rep( apply( ProbMat[ , idx ], 1, sum ), n )
                Mu[ , s ] <-  F1 / F2
            
                ## estimate sigma1
                F3 <- apply( Y * Y * crossprod( D, ProbMat[ , idx ] ), 1, sum )
                Vari <- F3 / F2 - Mu[ , s ] ^ 2
                Sigma[ , s] <- Mu[ , s ] / ( Vari / Mu[ , s ] - 1 )
                Sigma[ Sigma[ , s ] < 0 ] <- 100
            
            }
            
        } else {

            for( s in seq_len( S ) ){
                ## estimate for mu1
                idx <- seq_len( I ) + I * ( s - 1 )
                F1 <- apply( log( Y + 1 ) * crossprod( D, ProbMat[ , idx ] ), 1, sum )
                F2 <- apply( crossprod( D, ProbMat[ , idx ] ), 1, sum )
                ##rep( apply( ProbMat[ , idx ], 1, sum ), n )
                Mu[ , s ] <-  F1 / F2
            
                ## estimate sigma1
                F3 <- apply( log( Y + 1 ) * log( Y + 1 ) * crossprod( D, ProbMat[ , idx ] ), 1, sum )
                Sigma[ , s ] <- F3 / F2 - Mu[ , s ] ^ 2
            }
            
        }

        ## order the means
        od <-  apply( Mu, 1, order )
        Mu <- matrix( Mu[ cbind( rep( seq_len( N ), each = S ), c( od ) ) ], ncol = S, byrow = TRUE )
        Sigma <- matrix( Sigma[ cbind( rep( seq_len( N ), each = S ), c( od ) ) ], ncol = S, byrow = TRUE )
        
        ## convert everything to matrices
        B <- matrix( rep( b, each = K ), nrow = K )

        ## format ProbMat
        ProbMat.Format <- matrix( 0, nrow = K * S, ncol = I )
        for( s in 1:S )
          ProbMat.Format[ ( s - 1 ) * K + seq_len( K ) , ] <- ProbMat[ , seq_len( I ) + I * ( s - 1 ) ]
        ProbMat <- ProbMat.Format

        oldpar <- newpar
        newpar <- c( c( W ), probz, zeta, c( P ), c( Mu ), c( Sigma ) )
        
    }## finish outer loop

    conv <- FALSE
    if( outitr < maxitr )
        conv <- TRUE
    
  new( "MBASICFit",
      Theta = bestTheta,
      W = bestW,
      Z = predZ,
      b = bestb,
      aic = - 2 * tail( alllik, 1 ) + 2 * numpar,
      bic = - 2 * tail( alllik, 1 ) + log( N * I ) * numpar,
      aicc = -2 * tail( alllik, 1 ) + 2 * numpar + 2 * numpar * ( numpar + 1 ) / ( N * I - numpar - 1 ),
      alllik = alllik,
      lik = tail( alllik, 1 ),
      zeta = zeta,
      Mu = Mu,
      Sigma = Sigma,
      probz = probz,
      P = P,
      converged = conv,
      Theta.err = tail( allerr, 1 ),
      ARI = tail( allari, 1 ),
      W.err = tail( W.err, 1 ),
      MisClassRate = tail( allmisclass, 1 )
      )
}

