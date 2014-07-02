#include "mcmc.h"

SEXP mcmc( SEXP _b, SEXP _States, SEXP _Theta, SEXP _W, SEXP _P, SEXP _Mu, SEXP _Sigma, SEXP _D, SEXP _alpha, SEXP _betaw, SEXP _betap, SEXP _xi, SEXP _tau, SEXP _omega, SEXP _nu, SEXP _zeta, SEXP _Gamma, SEXP _Y) {

  // The following values are updated in MCMC iterations
  NumericVector b(_b); // length I
  NumericVector States(_States); // length I
  NumericMatrix Theta(_Theta); // I by KS
  NumericMatrix W(_W); // KS by I + 1
  NumericMatrix P(_P); // I by S
  NumericMatrix Mu(_Mu); // N by S
  NumericMatrix Sigma(_Sigma); // N by S
  NumericMatrix D(_D); // K by N
  double zeta = as<double>(_zeta); 

  // The following values are piror parameters and are fixed
  double alpha = as<double>(_alpha);
  double betaw = as<double>(_betaw);
  double betap = as<double>(_betap);
  double xi = as<double>(_xi);
  double tau = as<double>(_tau);
  double omega = as<double>(_omega);
  double nu = as<double>(_nu);

  // The following is the external information.
  NumericMatrix Gamma(_Gamma); // I by N
  NumericMatrix Y(_Y); // I by N

  // extract the dimensions
  int I = P.nrow();
  int S = P.ncol();
  int K = W.nrow();
  int N = D.ncol();

  double _LOW = 1e-10;

  // iterators
  int i, j, k, s;//, likid;

  // Compute J
  int J = 0;
  for(i = 0; i < I; i ++) {
    if(J < States[i]) {
      J = States[i];
    }
  }

  // Compute the density matrix
  NumericMatrix logF(I, K * S);

  // compute the density matrix
  for(i = 0; i < I; i ++) {
    for(k = 0; k < K; k ++) {
      for(s = 0; s < S; s ++) {
        logF(i, (s - 1) * K + k) = 0; // initialize
        for(n = 0; n < N; n ++) {
          if(D(k, n) == 1) {
            double tmp = R::dnorm(log(Y(i, n) + 1), Mu(n, s) * Gamma(i, n), Sigma(n, s), 0);
            if(tmp < _LOW) {
              tmp = _LOW;
            }
            logF(i, (s - 1) * K + k) = logF(i, (s - 1) * K + k) + log(tmp);
          }
        }
      }
    }
  }

  // Sample for b

  if(zeta < _LOW) {
    zeta = _LOW;
  } else if(zeta > 1 - _LOW) {
    zeta = 1 - _LOW;
  }

  for(i = 0; i < I; i ++) {
    // log probability of drawing 1
    double prob1 = 0;
    for(s = 0; s < S; s ++) {
      if(P(i, s) < _LOW)
        P(i, s) = _LOW;
      double tmpsum = 0;
      for(k = 0; k < K; k ++) {
        tmpsum += Theta(i, (s - 1) * K + k);
      }
      prob1 += tmpsum * log(P(i, s));
    }
    prob1 += log(zeta);

    // log probability of drawing 0

    double prob0 = 0;
    for(s = 0; s < S; s ++) {
      for(k = 0; k < K; k ++) {
        prob0 += Theta(i, (s - 1) * K + k) * log(W(k + (s - 1) * K, States[i]));
      }
    }
    prob0 += log(1 - zeta);

    // draw sample
    if(prob1 > prob0)
      b(i) = R::rbinom(1, 1 / (1 + exp(prob0 - prob1)));
    else
      b(i) = R::rbinom(1, exp(prob1 - prob0) / (1 + exp(prob1 - prob0)));
  }

  // Sample for States
  for(i = 0; i < I; i ++) {
    double probz[J + 1];
    // probability of choosing a existing sample
    for(j = 0; j < J; j ++) {
      double n_j = ClusterSize[j];// number of elements in this cluster
      if(States[i] == j)
        n_j --;
      if(n_j < _LOW) {
        n_j = _LOW; // approxmiation
      }
      probz[j] = log(n_j);
      for(k = 0; k < K; k ++) {
        for(s = 0; s < S; s ++) {
          probz[j] += Theta(i, (s - 1) * K + k) * (1 - b(i)) * log(W((s - 1) * K + k, j));
        }
      }
    }

    // create a new sample
    probz[J] = log(alpha);
    //notice: only calculate the Gamma functions if b_i = 0;
    if(b(i) == 0) {
      probz[J] -= log(S) * K;
    }

    double max_exp = probz[0];
    for(j = 1; j < J + 1; j ++) {
      if(max_exp < probz[j]) {
        max_exp = probz[j];
      }
    }

    double cum_probz[J + 1];
    cum_probz[0] = exp(probz[0] - maxexp);
    for(j = 1; j < J + 1; j ++ ) {
      cum_probz[j] = cumprobz[j - 1] + exp(probz[j] - maxexp);
    }

    double tmp = R::runif(0, cum_probz[J]);
    int selectj = J;
    while(selectj > 0 && tmp < cum_probz[selectj]) {
      selectj --;
    }
    if(cumprob[selectj] <= tmp)
      selectj ++;

    States[i] = selectj;

    // sample the new w
    NumericVector neww(K * S);
    for(k = 0; k < K; k ++) {
      double tmpGamma[S + 1];
      for(s = 0; s < S; s ++) {
        tmpGamma[s] = rgamma(betaw + (1 - b(i)) * Theta(i, (k - 1) * S + k), 1);
        tmpGamma[S] += tmpGamma[s];
      }
      for(s = 0; s < S; s ++) {
        neww[(s - 1) * K + k] = tmpGamma[s] / tmpGamma[S];
      }
    }

    W(_, J) = neww;

    //compute the number of units in each cluster
    int ClusterSize[J + 1];
    int NewLabel[J + 1];
    for(i = 0; i < I; i ++) {
      ClusterSize[States[i]] ++;
    }
    int newJ = 0;
    for(j = 0; j < J + 1; j ++) {
      if(ClusterSize[j] > 0) {
        NewLabel[j] = newJ;
        newJ ++;
      }
    }
    int OldLabel[newJ + 1];
    int newj = 0;
    for(j = 0; j < J + 1; j ++) {
      if(ClusterSize[j] > 0) {
        OldLabel[newj] = j;
        newj ++;
      }
    }

    //Use new labels for each unit
    for(i = 0; i < I; i ++) {
      States[i] = NewLabel[States[i]];
    }

    //Clean columns in W
    for(j = 0; j < newJ; j ++) {
      W(_, j) = W(_, OldLabel[j]);
    }
    for(k = 0; k < K; k ++)
      for(s = 0; s < S; s ++)
        W((s - 1) * K + k, newJ + 1) = 0;

  }

  Rcpp::List ret = Rcpp::List::create(
      Rcpp::Named("P") = Theta_b,
      Rcpp::Named("zeta") = zeta,
      Rcpp::Named("probz") = Z_mean,
      Rcpp::Named("W") = W_max,
      Rcpp::Named("Theta_mean") = Theta_mean,
      Rcpp::Named("b_prob") = b_mean,
      Rcpp::Named("predZ") = predZ
      //,Rcpp::Named("oldlik") = oldlik

  return( ret );

}
