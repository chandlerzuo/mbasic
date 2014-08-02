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
  NumericMatrix W(K * S, I + 1);
  NumericMatrix P(I, S);

  // iterators
  int i, j, k = 0, s, n, i1;//, likid;
  double loss;

  double _LOW = 1e-10;

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
      double tmp[S + 1];
      tmp[S] = 0;
      for(s = 0; s < S; s ++) {
	tmp[s] = _LOW;
	tmp[S] += tmp[s];
	W(s * K + k, j) = 0;
      }
      for(i = 0; i < I; i ++) {
	if(b[i] == 0 && States[i] == j) {
	  tmp[Theta(i, k)] ++;
	  tmp[S] ++;
	}
      }
      for(s = 1; s < S; s ++) {
	      W(s * K + k, j) = tmp[s] / tmp[S];
      }
    }
  }

  // update P
  for(i = 0; i < I; i ++) {
    double tmp[S + 1];
    tmp[S] = 0;
    for(s = 0; s < S; s++) {
      tmp[s] = _LOW;
      tmp[S] += _LOW;
      P(i, s) = 0;
    }
    for(k = 0; k < K; k ++) {
      tmp[Theta(i, k)] ++;
      tmp[S] ++;
    }
    for(s = 0; s < S; s ++) {
	    P(i, s) = tmp[s] / tmp[S];
    }
  }

  // update b
  for(i = 0; i < I; i ++) {
    double tmp = 0;
    for(k = 0; k < K; k ++) {
	    tmp += 2 * lambdap * (1 - P(i, Theta(i, k)));
    }

    for(k = 0; k < K; k ++) {
	    tmp -= 2 * lambdaw * (1 - W(Theta(i, k) * K + k, States[i]);
    }

    if(tmp < 0)
      b[i] = 1;
    else
      b[i] = 0;
  }

  loss = 0;
  for(i = 0; i < I; i ++) {
    for(n = 0; n < N; n ++) {
      k = D[n];
      s = Theta(i, k);
      loss += (Y(i, n) - Mu(n, s) * Gamma(i, s * N + n)) *
	(Y(i, n) - Mu(n, s) * Gamma(i, s * N + n));
    }

    if(b[i] == 1) {
      for(k = 0; k < K; k ++) {
	      loss += 2 * lambdap * (1 - P(i, Theta(i, k)));
      }
    } else {
      for(k = 0; k < K; k ++) {
	      loss += 2 * (1 - W(Theta(i, k) * K + k, States[i])) * lambdaw;
      }
    }
  }
  loss += lambda * (J - 1);

  printf("b, Loss function = %3.3f, number of clusters = %d\n", loss, J);

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
            tmp[s] += (Y(i, n) - Mu(n, s) * Gamma(i, s * N + n)) * (Y(i, n) - Mu(n, s) * Gamma(i, s * N + n));
          }
        }
        if(b[i] == 1) {
		tmp[s] += 2 * lambdap * (1 - P(i, Theta(i, k)));
	} else {
		tmp[s] += 2 * lambdaw * (1 - W(Theta(i, k) * K + k, States[i]));
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

  loss = 0;
  for(i = 0; i < I; i ++) {
    for(n = 0; n < N; n ++) {
      k = D[n];
      s = Theta(i, k);
      loss += (Y(i, n) - Mu(n, s) * Gamma(i, s * N + n)) *
	(Y(i, n) - Mu(n, s) * Gamma(i, s * N + n));
    }

    if(b[i] == 1) {
      for(k = 0; k < K; k ++) {
	      loss += 2 * lambdap * (1 - P(i, Theta(i, k)));
      }
    } else {
      for(k = 0; k < K; k ++) {
	      loss += 2 * (1 - W(Theta(i, k) * K + k, States[i])) * lambdaw;
      }
    }
  }
  loss += lambda * (J - 1);

  printf("theta, Loss function = %3.3f, number of clusters = %d\n", loss, J);

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
      for(k = 0; k < K; k ++) {
	      tmp[j] += 2 * (1 - W(Theta(i, k) * K + k, States[i])) * lambdaw;
      }
      // assign the minimum cost
      if(tmp[j] < mintmp) {
        mintmp = tmp[j];
        newState = j;
      }
    }

    States[i] = newState;
    ClusterSize[newState] ++;
    if(mintmp >= lambda) {
      // a new cluster is formed
      if(J != newState) {
        printf("Error: new state is not J = %d.", J);
        exit(1);
      }
      J ++;
      for(s = 0; s < S; s ++) {
	      for(k = 0; k < K; k ++) {
		      if(Theta(i, k) != s) {
			      W(s * K + k, J) = _LOW;
		      } else {
			      W(s * K + k, J) = 1 - (S - 1) * _LOW;
		      }
	      }
      }
    }
    if(ClusterSize[oldState] == 0) {
      // an old cluster should be removed
      W(_, oldState) = W(_, J - 1);
      for(k = 0; k < K * S; k ++) {
        W(k, J - 1) = 0;
      }
      for(i1 = 0; i1 < I; i1 ++) {
        if(States[i1] == J - 1)
          States[i1] = oldState;
      }
      ClusterSize[oldState] = ClusterSize[J - 1];
      ClusterSize[J - 1] = 0;
      J --;
    }
  }

  loss = 0;
  for(i = 0; i < I; i ++) {
    for(n = 0; n < N; n ++) {
      k = D[n];
      s = Theta(i, k);
      loss += (Y(i, n) - Mu(n, s) * Gamma(i, s * N + n)) *
	(Y(i, n) - Mu(n, s) * Gamma(i, s * N + n));
    }

    if(b[i] == 1) {
      for(k = 0; k < K; k ++) {
	      loss += 2 * lambdap * (1 - P(i, Theta(i, k)));
      }
    } else {
      for(k = 0; k < K; k ++) {
	      loss += 2 * (1 - W(Theta(i, k) * K + k, States[i])) * lambdaw;
      }
    }
  }
  loss += lambda * (J - 1);

  printf("z, Loss function = %3.3f, number of clusters = %d\n", loss, J);

  // update mu
  for(n = 0; n < N; n ++) {
    for(s = 0; s < S; s ++) {
      double denom = 0, numer = 0;
      for(i = 0; i < I; i ++) {
	if(Theta(i, D[n]) == s) {
	  numer += Y(i, n);
	  denom += Gamma(i, n + s * N);
	}
      }
      Mu(n, s) = numer / denom;
    }
  }

  // calculate the loss function

  loss = 0;
  for(i = 0; i < I; i ++) {
    for(n = 0; n < N; n ++) {
      k = D[n];
      s = Theta(i, k);
      loss += (Y(i, n) - Mu(n, s) * Gamma(i, s * N + n)) *
	(Y(i, n) - Mu(n, s) * Gamma(i, s * N + n));
    }

    if(b[i] == 1) {
      for(k = 0; k < K; k ++) {
	      loss += 2 * lambdap * (1 - P(i, Theta(i, k)));
      }
    } else {
      for(k = 0; k < K; k ++) {
	      loss += 2 * (1 - W(Theta(i, k) * K + k, States[i])) * lambdaw;
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
