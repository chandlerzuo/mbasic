#include "mcmc.h"

//this file is in development

SEXP mcmc( SEXP _b, SEXP _States, SEXP _Theta, SEXP _Mu, SEXP _Sigma, SEXP _D, SEXP _alpha, SEXP _betaw, SEXP _betap, SEXP _xi, SEXP _tau, SEXP _omega, SEXP _nu, SEXP _zeta, SEXP _Gamma, SEXP _Y) {

	// The following values are updated in MCMC iterations
	NumericVector b(_b); // length I
	IntegerVector States(_States); // length I
	IntegerMatrix Theta(_Theta); // I by K
	//NumericMatrix W(_W); // KS by I + 1
	//NumericMatrix P(_P); // I by S
	NumericMatrix Mu(_Mu); // N by S
	NumericMatrix Sigma(_Sigma); // N by S
	IntegerVector D(_D); // Length N, valued in {0, 1, ..., K-1}
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

	double _LOW = 1e-10;

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

	// Compute the density matrix
	NumericMatrix logF(I, K * S);

	// compute the density matrix
	for(i = 0; i < I; i ++) {
		for(k = 0; k < K; k ++) {
			for(s = 0; s < S; s ++) {
				logF(i, s * K + k) = 0; // initialize
				for(n = 0; n < N; n ++) {
					if(D[n] == k) {
						double tmp = R::dnorm(log(Y(i, n) + 1), Mu(n, s) * Gamma(i, s * N + n), Sigma(n, s), 0);
						if(tmp < _LOW) {
							tmp = _LOW;
						}
						logF(i, s * K + k) += log(tmp);
					}
				}
			}
		}
	}
	printf("Finished computing the density matrix\n");

	// Update W
	for(j = 0; j < J; j ++) {
		for(k = 0; k < K; k ++) {
			double betaRates[S], randGamma[S + 1];
			for(i = 0; i < I; i ++) {
				if(j == States[i] && b[i] == 0) {
					betaRates[Theta(i, k)] ++;
				}
			}
			for(s = 0; s < S; s ++) {
				betaRates[s] += betaw;
				randGamma[s] = R::rgamma(betaRates[s], 1);
				randGamma[S] += randGamma[s];
			}
			for(s = 0; s < S; s ++) {
				W(s * K + k, j) = randGamma[s] / randGamma[S];
				if(W(s * K + k, j) < _LOW) {
				  W(s * K + k, j) = _LOW;
				}
			}
		}
	}
	printf("Finished updating W\n");
	
	// update for P
	for(i = 0; i < I; i ++) {
		double betaRates[S], randGamma[S + 1];
		for(k = 0; k < K; k ++) {
			if(b[i] == 1) {
				betaRates[Theta(i, k)] ++;
			}
		}
		for(s = 0; s < S; s ++) {
			betaRates[s] += betap;
			randGamma[s] = R::rgamma(betaRates[s], 1);
			randGamma[S] += randGamma[s];
		}
		for(s = 0; s < S; s ++) {
			P(i, s) = randGamma[s] / randGamma[S];
			if(P(i, s) < _LOW) {
			  P(i, s) = _LOW;
			}
		}
	}
	printf("Finished updating P\n");

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
                                 if(Theta(i, k) == s) {
                                         tmpsum ++;
                                 }
                         }
                         prob1 += tmpsum * log(P(i, s));
                 }
                 prob1 += log(zeta);

                 // log probability of drawing 0

                 double prob0 = 0;
                 for(k = 0; k < K; k ++) {
                         prob0 += log(W(k + (Theta(i, k) - 1) * K, States[i]));
                 }

                 prob0 += log(1 - zeta);

                 // draw sample
                 if(prob1 > prob0)
                         b(i) = R::rbinom(1, 1 / (1 + exp(prob0 - prob1)));
                 else
                         b(i) = R::rbinom(1, exp(prob1 - prob0) / (1 + exp(prob1 - prob0)));
         }
         printf("Finished updating b\n");

         // Sample for States
         for(i = 0; i < I; i ++) {
                 double probz[J + 1];
                 // probability of choosing a existing sample
                 ClusterSize[States[i]] --;
                 int oldState = States[i];
                 for(j = 0; j < J; j ++) {
                         double n_j = ClusterSize[j];// number of elements in this cluster
                         if(n_j < _LOW) {
                                 n_j = _LOW; // approxmiation
                         }
                         probz[j] = log(n_j);
                         for(k = 0; k < K; k ++) {
                                 probz[j] += (1 - b(i)) * log(W((Theta(i, k) - 1) * K + k, j));
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

                 double cum_probz[J + 2];
                 cum_probz[0] = 0;
                 for(j = 0; j < J + 1; j ++ ) {
                         cum_probz[j + 1] = cum_probz[j] + exp(probz[j] - max_exp);
                 }

                 double tmp = R::runif(0, cum_probz[J + 1]);
                 int selectj = J + 1;
                 while(tmp <= cum_probz[selectj]) {
                         selectj --;
                 }

                 States[i] = selectj;
                 ClusterSize[selectj] ++;

                 if(selectj == J) {
                         // sample the new w
                         NumericVector neww(K * S);
                         for(k = 0; k < K; k ++) {
                                 double tmpGamma[S + 1];
                                 for(s = 0; s < S; s ++) {
                                         if(s == Theta(i, k) && b[i] == 0) {
                                           tmpGamma[s] = R::rgamma(betaw + 1, 1);
                                         } else {
                                           tmpGamma[s] = R::rgamma(betaw, 1);
                                         }
                                         tmpGamma[S] += tmpGamma[s];
                                 }
                                 for(s = 0; s < S; s ++) {
                                         neww[s * K + k] = tmpGamma[s] / tmpGamma[S];
                                 }
                         }

                         W(_, J) = neww;
                 }
                 //		printf("Old State = %d, new state = %d, old state cluster size = %d, new cluster size = %d.\n", oldState, selectj, ClusterSize[oldState], ClusterSize[selectj]);

                 // switch cluster labels if the cluster size changes
                 if(ClusterSize[oldState] == 0) {
                         if(selectj == J) {
                                 if(ClusterSize[selectj] != 1) 
                                         printf("Error: new cluster size is not 1, J = %d, selectj = %d, oldState = %d\n", J, selectj, oldState);
                                 States[i] = oldState;
                                 W(_, oldState)  = W(_, selectj);
                                 ClusterSize[oldState] = ClusterSize[selectj];
                                 ClusterSize[selectj] = 0;
                         } else {
                                 for(i = 0; i < I; i ++) {
                                         if(States[i] == J - 1) {
                                                 States[i] = oldState;
                                                 if(States[i] > J - 1) {
                                                   printf("Error: some cluster label are larger than J - 1, i = %d, State = %d\n", i, States[i]);
                                                   return(wrap(-1));
                                                 }
                                         }
                                 }
                                 W(_, oldState) = W(_, J - 1);
                                 for(i = 0; i < K * S; i ++) {
                                         W(i, J - 1) = 0;
                                 }
                                 ClusterSize[oldState] = ClusterSize[J - 1];
                                 ClusterSize[J - 1] = 0;
                                 J--;
                         }
                 } else if(selectj == J) {
                         J ++;
                 }

         }
         printf("Finished updating States\n");

         // update for Theta
         for(i = 0; i < I; i ++) {
                 for(k = 0; k < K; k ++) {
                         double logProb[S];
                         for(s = 0; s < S; s ++) {
                                 logProb[s] = logF(i, s * K + k);
                                 if(b[i] == 0) {
                                         logProb[s] += log(W(s * K + k, States[i]));
                                 } else {
					logProb[s] += log(P(i, s));
			  }
			}
			double maxLogProb = logProb[0];
			for(s = 1; s < S; s ++) {
				if(maxLogProb < logProb[s])
					maxLogProb = logProb[s];
			}
			double cumProb[S + 1];
			cumProb[0] = 0;
			for(s = 0; s < S; s ++) {
				cumProb[s + 1] = cumProb[s] + exp(logProb[s] - maxLogProb);
			}
			double rndVar = R::runif(0, cumProb[S]);
			s = S;
			while(rndVar <= cumProb[s]) {
				s--;
			}
			Theta(i, k) = s;
			if(s == S) {
				printf("Error: sampling state can not be S, i = %d, k = %d, rndVar = %3.3f.\n", i, k, rndVar);
				for(s = 0; s < S + 1; s ++) {
					printf("Cumprob[s] = %3.3f\t", cumProb[s]);
				}
				printf("\n");
			}
		}
	}

	printf("Finished updating Theta\n");

	// update for mu and sigma
	for(n = 0; n < N; n ++) {
		for(s = 0; s < S; s ++) {
			double productTerm = 0, squareTerm = 0, varTerm = 0;
			for(i = 0; i < I; i ++) {
				productTerm += Gamma(i, s * K + k) * log(Y(i, n) + 1);
				squareTerm += Gamma(i, s * K + k) * Gamma(i, s * K + k);
				varTerm += (log(Y(i, n) + 1) - Mu(n, s) * Gamma(i, s * K + k)) * (log(Y(i, n) + 1) - Mu(n, s) * Gamma(i, s * K + k)) ;
			}
			Mu(n, s) = R::rnorm((tau * productTerm + xi * Sigma(n, s)) / (tau * squareTerm + Sigma(n, s)), Sigma(n, s) * tau / (tau * squareTerm + Sigma(n, s)));
			Sigma(n, s) = 1 / R::rgamma(omega + 1 + I / 2,
						 nu + varTerm / 2);
		}
	}
	
	printf("Finished updating mu and sigma\n");

	// update zeta
	zeta = 0;
	for(i = 0; i < I; i ++) {
		zeta += b[i];
	}
	zeta /= I;

	Rcpp::List ret = Rcpp::List::create(
					    Rcpp::Named("Theta") = Theta,
					    Rcpp::Named("States") = States,
					    Rcpp::Named("b") = b,
					    Rcpp::Named("W") = W,
					    Rcpp::Named("P") = P,
					    Rcpp::Named("Mu") = Mu,
					    Rcpp::Named("Sigma") = Sigma,
					    Rcpp::Named("zeta") = zeta
					    //,Rcpp::Named("oldlik") = oldlik
					    );
	return( ret );
  
}
