#include "mcmc.h"

//this file is in development

SEXP map( SEXP _b, SEXP _States, SEXP _Theta, SEXP _Mu, SEXP _D,
	  SEXP _Gamma, SEXP _Y, SEXP _lambdap, SEXP _lambdaw, SEXP _lambda) {

  // The following values are 1updated in MCMC iterations
  IntegerVector b(_b); // length I
  IntegerVector States(_States); // length I
  IntegerMatrix Theta(_Theta); // I by K
  //NumericMatrix W(_W); // KS by I + 1
  //NumericMatrix P(_P); // I by S
  NumericMatrix Mu(_Mu); // N by S
  IntegerVector D(_D); // Length N, valued in {0, 1, ..., K-1}
  double lambda = as<double>(_lambda);

  // The following values are piror parameters and are fixed
  double lambdap = as<double>(_lambdap);
  double lambdaw = as<double>(_lambdaw);

  // The following is the external information.
  NumericMatrix Gamma(_Gamma); // I by N * S
  NumericMatrix Y(_Y); // I by N

  // extract the dimensions
  int I = b.size();
  int S = Mu.ncol();
  int K = Theta.ncol();
  int N = D.size();

  // The following will be computed
  IntegerMatrix W(K, I + 1);
  IntegerVector P(I);

  // iterators
  int i, j, k = 0, s, n;//, likid;

  int ClusterSize[I + 1];

  for(i = 0; i < I + 1; i ++) {
    ClusterSize[i] = 0;
  }

  // Compute J

  int J = 0;
  for(i = 0; i < I; i ++) {
    if(J < States[i]) {
      J = States[i];
    }
    ClusterSize[States[i]] ++;
  }
  J ++;

  // Update W
  for(j = 0; j < J; j ++) {
    for(k = 0; k < K; k ++) {
      double tmp[S];
      for(s = 0; s < S; s ++) {
	tmp[s] = 0;
      }
      for(i = 0; i < I; i ++) {
	if(b[i] == 0 && States[i] == j) {
	  tmp[Theta(i, k)] ++;
	}
      }
      W(k, j) = 0;
      for(s = 1; s < S; s ++) {
	if(tmp[W(k, j)] > tmp[s]) {
	  W(k, j) = s;
	}
      }
    }
  }

  // update P
  for(i = 0; i < I; i ++) {
    double tmp[S];
    for(s = 0; s < S; s++) {
      tmp[s] = 0;
    }
    for(k = 0; k < K; k ++) {
      tmp[Theta(i, k)] ++;
    }
    P[i] = 0;
    for(s = 1; s < S; s ++) {
      if(tmp[P[i]] < tmp[s]) {
        P[i] = s;
      }
    }
  }

  // update b
  for(i = 0; i < I; i ++) {
    double tmp = 0;
    for(k = 0; k < K; k ++) {
      if(Theta(i, k) != P[i]) {
        tmp += lambdap;
      }
    }

    for(j = 0; j < J; j ++) {
      if(States[i] == j) {
        for(k = 0; k < K; k ++) {
          if(W(k, j) != Theta(i, k)) {
            tmp -= lambdaw;
          }
        }
      }
    }

    if(tmp < 0)
      b[i] = 1;
    else
      b[i] = 0;
  }

  // update Theta
  for(i = 0; i < I; i ++) {
    for(k = 0; k < K; k ++) {
      double tmp[S];
      for(s = 0; s < S; s ++) {
        // initialize
        tmp[s] = 0;
        //
        for(n = 0; n < N; n ++) {
          if(D[n] == k) {
            tmp[s] += (log(Y(i, n) + 1) - Mu(n, s) * Gamma(i, s * N + n)) * (log(Y(i, n) + 1) - Mu(n, s) * Gamma(i, s * N + n));
          }
        }
        if(s != P[i] && b[i] == 1)
          tmp[s] += lambdap;
        if(b[i] == 0 && s != W(k, States[i])) {
          tmp[s] += lambdaw;
        }
      }
      // Assign new values
      Theta(i, k) = 0;
      for(s = 1; s < S; s ++) {
        if(tmp[s] < tmp[Theta(i, k)])
          Theta(i, k) = s;
      }
    }
  }


  //update clusters
  for(i = 0; i < I; i ++) {
    // do not update cluster if this is a singleton
    if(b[i] == 1)
      continue;

    int oldState = States[i];
    ClusterSize[States[i]] --;
    double tmp[J];
    double mintmp = lambda;// cost of starting a new cluster
    int newState = J;
    for(j = 0; j < J; j ++) {
      tmp[j] = 0;
      if(b[i] == 0) {
        for(k = 0; k < K; k ++) {
          if(Theta(i, k) != W(k, j)) {
            tmp[j] ++;
          }
        }
      }
      // assign the minimum cost
      if(tmp[j] < mintmp) {
        mintmp = tmp[j];
        newState = j;
      }
    }

    States[i] = newState;
    ClusterSize[newState] ++;
    if(mintmp > lambda) {
      // a new cluster is formed
      if(J != newState) {
        printf("Error: new state is not J = %d.", J);
        exit(1);
      }
      J ++;
      W(_, newState) = Theta(i, _);
    }
    if(ClusterSize[oldState] == 0) {
      // an old cluster should be removed
      W(_, oldState) = W(_, J - 1);
      for(k = 0; k < K; k ++) {
        W(k, J - 1) = 0;
      }
      for(i = 0; i < I; i ++) {
        if(States[i] == J - 1)
          States[i] = oldState;
      }
      ClusterSize[oldState] = ClusterSize[J - 1];
      ClusterSize[J - 1] = 0;
      J --;
    }
  }

  // update mu
  for(n = 0; n < N; n ++) {
    for(s = 0; s < S; s ++) {
      double denom = 0, numer = 0;
      for(i = 0; i < I; i ++) {
        numer += log(Y(i, n) + 1);
        denom += Gamma(i, n + s * N);
      }
      Mu(n, s) = numer / denom;
    }
  }

  // calculate the loss function

  double loss = 0;
  for(i = 0; i < I; i ++) {
    for(n = 0; n < N; n ++) {
      k = D[n];
      s = Theta(i, k);
      loss += (log(Y(i, n) + 1) - Mu(n, s) * Gamma(i, s * N + n)) *
	(log(Y(i, n) + 1) - Mu(n, s) * Gamma(i, s * N + n));
    }
    
    if(b[i] == 1) {
      for(k = 0; k < K; k ++) {
	if(Theta(i, k) != P[i]) {
	  loss += lambdap;
	}
      }
    } else {
      for(k = 0; k < K; k ++) {
	if(W(k, States[i]) != Theta(i, k)) {
	  loss += lambdaw;
	}
      }
    }
  }
  loss += lambda * (J - 1);
  
  printf("Loss function = %3.3f, number of clusters = %d\n", loss, J);

  Rcpp::List ret = Rcpp::List::create(
				      Rcpp::Named("Theta") = Theta,
				      Rcpp::Named("States") = States,
				      Rcpp::Named("b") = b,
				      Rcpp::Named("W") = W,
				      Rcpp::Named("P") = P,
				      Rcpp::Named("Mu") = Mu,
				      Rcpp::Named("loss") = loss
				      //,Rcpp::Named("oldlik") = oldlik
				      );
  return( ret );

}
