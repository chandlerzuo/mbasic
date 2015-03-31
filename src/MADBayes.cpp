#include "mcmc.h"

//this file is under testing

SEXP MADBayes( SEXP _b, SEXP _clusterLabels, SEXP _Theta, SEXP _Mu, SEXP _D,
	       SEXP _Gamma, SEXP _Y, SEXP _lambdap, SEXP _lambdaw, SEXP _lambda) {

	// The following values are 1updated in MCMC iterations
	IntegerVector b(_b); // length I
	IntegerVector clusterLabels(_clusterLabels); // length I
	IntegerMatrix Theta(_Theta); // K by I
	//NumericMatrix W(_W); // KS by I + 1
	//NumericMatrix P(_P); // I by S
	NumericMatrix Mu(_Mu); // N by S
	IntegerVector D(_D); // Length N, valued in {0, 1, ..., K-1}
	double lambda = as<double>(_lambda);

	// The following values are piror parameters and are fixed
	double lambdap = as<double>(_lambdap);
	double lambdaw = as<double>(_lambdaw);

	// The following is the external information.
	NumericMatrix Gamma(_Gamma); // N*S by I
	NumericMatrix Y(_Y); // N by I

	// extract the dimensions
	int I = b.size();
	int S = Mu.ncol();
	int K = Theta.nrow();
	int N = D.size();

	// The following will be computed
	NumericMatrix W(K * S, I + 1);
	NumericMatrix P(I, S);

	// iterators
	int i, j, k = 0, s = 0, n, i1;//, likid;
	double loss;

	double _LOW = 1e-10;

	int ClusterSize[I + 1];

	for(i = 0; i < I + 1; i ++) {
		ClusterSize[i] = 0;
	}

	// Compute J

	int J = 0;
	for(i = 0; i < I; i ++) {
		if(J < clusterLabels[i]) {
			J = clusterLabels[i];
		}
		ClusterSize[clusterLabels[i]] ++;
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
				if(b[i] == 0 && clusterLabels[i] == j) {
					tmp[Theta(k, i)] ++;
					tmp[S] ++;
				}
			}
			for(s = 0; s < S; s ++) {
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
			tmp[Theta(k, i)] ++;
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
			tmp += 2 * lambdap * (1 - P(i, Theta(k, i)));
		}

		for(k = 0; k < K; k ++) {
			tmp -= 2 * lambdaw * (1 - W(Theta(k, i) * K + k, clusterLabels[i]));
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
			s = Theta(k, i);
			loss += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) *
				(Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
		}

		if(b[i] == 1) {
			for(k = 0; k < K; k ++) {
				loss += 2 * lambdap * (1 - P(i, Theta(k, i)));
			}
		} else {
			for(k = 0; k < K; k ++) {
				loss += 2 * (1 - W(Theta(k, i) * K + k, clusterLabels[i])) * lambdaw;
			}
		}
	}
	loss += lambda * (J - 1);

	//printf("b, Loss function = %3.3f, number of clusters = %d\n", loss, J);

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
						tmp[s] += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) * (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
					}
				}
				if(b[i] == 1) {
					tmp[s] += 2 * lambdap * (1 - P(i, s));
				} else {
					tmp[s] += 2 * lambdaw * (1 - W(s * K + k, clusterLabels[i]));
				}
			}
			// Assign new values
			Theta(k, i) = 0;
			for(s = 1; s < S; s ++) {
				if(tmp[s] < tmp[Theta(k, i)])
					Theta(k, i) = s;
			}
		}
	}

	loss = 0;
	for(i = 0; i < I; i ++) {
		for(n = 0; n < N; n ++) {
			k = D[n];
			s = Theta(k, i);
			loss += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) *
				(Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
		}

		if(b[i] == 1) {
			for(k = 0; k < K; k ++) {
				loss += 2 * lambdap * (1 - P(i, Theta(k, i)));
			}
		} else {
			for(k = 0; k < K; k ++) {
				loss += 2 * (1 - W(Theta(k, i) * K + k, clusterLabels[i])) * lambdaw;
			}
		}
	}
	loss += lambda * (J - 1);

	//printf("theta, Loss function = %3.3f, number of clusters = %d\n", loss, J);

	//update clusters
	for(i = 0; i < I; i ++) {
		// do not update cluster if this is a singleton
		if(b[i] == 1)
			continue;

		int oldState = clusterLabels[i];
		ClusterSize[clusterLabels[i]] --;
		double tmp[J];
		double mintmp = lambda;// cost of starting a new cluster
		int newState = J;
		for(j = 0; j < J; j ++) {
			tmp[j] = 0;
			for(k = 0; k < K; k ++) {
				tmp[j] += 2 * (1 - W(Theta(k, i) * K + k, j)) * lambdaw;
			}
			// assign the minimum cost
			if(tmp[j] < mintmp) {
				mintmp = tmp[j];
				newState = j;
			}
		}

		clusterLabels[i] = newState;
		ClusterSize[newState] ++;
		if(mintmp >= lambda) {
			// a new cluster is formed
			if(J != newState) {
				printf("Error: new state is not J = %d.", J);
				exit(1);
			}
			for(s = 0; s < S; s ++) {
				for(k = 0; k < K; k ++) {
					if(Theta(k, i) != s) {
						W(s * K + k, J) = _LOW;
					} else {
						W(s * K + k, J) = 1 - (S - 1) * _LOW;
					}
				}
			}
			J ++;
		}
		if(ClusterSize[oldState] == 0) {
			// an old cluster should be removed
			W(_, oldState) = W(_, J - 1);
			for(k = 0; k < K * S; k ++) {
				W(k, J - 1) = 0;
			}
			for(i1 = 0; i1 < I; i1 ++) {
				if(clusterLabels[i1] == J - 1)
					clusterLabels[i1] = oldState;
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
			s = Theta(k, i);
			loss += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) *
				(Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
		}

		if(b[i] == 1) {
			for(k = 0; k < K; k ++) {
				loss += 2 * lambdap * (1 - P(i, Theta(k, i)));
			}
		} else {
			for(k = 0; k < K; k ++) {
				loss += 2 * (1 - W(Theta(k, i) * K + k, clusterLabels[i])) * lambdaw;
			}
		}
	}
	loss += lambda * (J - 1);

	//printf("z, Loss function = %3.3f, number of clusters = %d\n", loss, J);

	// update mu
	for(n = 0; n < N; n ++) {
		double denom = 0, numer = 0;
		for(i = 0; i < I; i ++) {
			denom += Gamma(n + s * N, i);
			numer += Y(n, i);
		}
		for(s = 0; s < S; s ++) {
			Mu(n, s) = numer / denom;
		}
		for(s = 0; s < S; s ++) {
			denom = 0;
			numer = 0;
			for(i = 0; i < I; i ++) {
				if(Theta(D[n], i) == s) {
					numer += Y(n, i);
					denom += Gamma(n + s * N, i);
				}
			}
			if(denom > 0) {
				Mu(n, s) = numer / denom;
			}
		}
	}

	// calculate the loss function

	loss = 0;
	for(i = 0; i < I; i ++) {
		for(n = 0; n < N; n ++) {
			k = D[n];
			s = Theta(k, i);
			loss += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) *
				(Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
		}

		if(b[i] == 1) {
			for(k = 0; k < K; k ++) {
				loss += 2 * lambdap * (1 - P(i, Theta(k, i)));
			}
		} else {
			for(k = 0; k < K; k ++) {
				loss += 2 * (1 - W(Theta(k, i) * K + k, clusterLabels[i])) * lambdaw;
			}
		}
	}
	loss += lambda * (J - 1);

	//printf("Loss function = %3.3f, number of clusters = %d\n", loss, J);

	Rcpp::List ret = Rcpp::List::create(
					    Rcpp::Named("Theta") = Theta,
					    Rcpp::Named("clusterLabels") = clusterLabels,
					    Rcpp::Named("b") = b,
					    Rcpp::Named("W") = W,
					    Rcpp::Named("P") = P,
					    Rcpp::Named("Mu") = Mu,
					    Rcpp::Named("loss") = loss
					    //,Rcpp::Named("oldlik") = oldlik
					    );
	return( ret );

}
