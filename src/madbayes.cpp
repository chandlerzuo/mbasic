#include "mcmc.h"

//this file is under testing

SEXP madbayes(SEXP _clusterLabels, SEXP _Theta, SEXP _Mu, SEXP _D,
	      SEXP _Gamma, SEXP _Y, SEXP _lambdaw, SEXP _lambda,
	      SEXP _maxitr, SEXP _tol) {
	
	// The following values are 1updated in MCMC iterations
	IntegerVector clusterLabels(_clusterLabels); // length I
	IntegerMatrix Theta(_Theta); // K by I
	NumericMatrix Mu(_Mu); // N by S
	IntegerVector D(_D); // Length N, valued in {0, 1, ..., K-1}
	double lambda = as<double>(_lambda);

	// The following values are piror parameters and are fixed
	double lambdaw = as<double>(_lambdaw);

	double tol = as<double>(_tol);
	int maxitr = as<int>(_maxitr);

	// The following is the external information.
	NumericMatrix Gamma(_Gamma); // N*S by I
	NumericMatrix Y(_Y); // N by I

	// extract the dimensions
	int I = Theta.ncol();
	int S = Mu.ncol();
	int K = Theta.nrow();
	int N = D.size();

	// The following will be computed
	NumericMatrix W(K * S, I + 1);
	NumericMatrix P(I, S);
	NumericMatrix Sigma(N, S);

	// Distance from each unit to each cluster
	NumericMatrix Dist(I, I + 1);
	IntegerVector clusterSizes(I + 1);

	// iterators
	int i, j, k = 0, s = 0, n, i1, j1, itr;//, likid;
	int firstLabel, lastLabel;
	double loss = 0, oldloss;

	double _LOW = 1e-10;

	for(i = 0; i < I + 1; i ++) {
		clusterSizes[i] = 0;
	}

	// Compute J
	int J = 0;
	int nextClusterLabel = 0;
	for(i = 0; i < I; i ++) {
		if(J < clusterLabels[i]) {
			J = clusterLabels[i];
		}
		clusterSizes[clusterLabels[i]] ++;
	}
	J ++;
	nextClusterLabel = J;
	// printf("start for lambda %3.3f\n", lambda / lambdaw);

	for(itr = 0; itr < maxitr; itr ++) {
		oldloss = loss;
		// Update W
		NumericVector counts(S);
		NumericVector w_sol(S);
		for(j = 0; j < I + 1; j ++) {
			if(clusterSizes[j] > 0) {
				for(k = 0; k < K; k ++) {
					for(s = 0; s < S; s ++) {
						counts[s] = 0;
					}
					for(i = 0; i < I; i ++) {
						if(clusterLabels[i] == j) {
							for(s = 0; s < S; s ++) {
								if(Theta(k, i) == s) {
									counts[s] ++;
								}
							}
						}
					}
					for(s = 0; s < S; s ++) {
						w_sol[s] = (counts[s] + _LOW) / ((double) clusterSizes[j] + S * _LOW);
					}
					for(s = 0; s < S; s ++) {
						W(s * K + k, j) = w_sol(s);
					}
				}
			}
		}

		//loss = ComputeLoss(D, Theta, Y, Mu, Gamma, W, clusterLabels, clusterSizes, lambda, lambdaw);
		//printf("update W, Loss function = %3.3f, number of clusters = %d\n", loss, J);

		// update Theta
		for(i = 0; i < I; i ++) {
			for(k = 0; k < K; k ++) {
				double tmp[S];
				for(s = 0; s < S; s ++) {
					// initialize
					tmp[s] = 0;
					for(n = 0; n < N; n ++) {
						if(D[n] == k) {
							tmp[s] += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) * (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
						}
					}
					tmp[s] += lambdaw * ((1 - W(s * K + k, clusterLabels[i])) * (1 - W(s * K + k, clusterLabels[i])) - W(s * K + k, clusterLabels[i]) * W(s * K + k, clusterLabels[i]));
				}
				// Assign new values
				Theta(k, i) = 0;
				for(s = 1; s < S; s ++) {
					if(tmp[s] < tmp[Theta(k, i)])
						Theta(k, i) = s;
				}
			}
		}

		//loss = ComputeLoss(D, Theta, Y, Mu, Gamma, W, clusterLabels, clusterSizes, lambda, lambdaw);
		//printf("update states, Loss function = %3.3f, number of clusters = %d\n", loss, J);

		//update clusters
		for(i = 0; i < I; i ++) {
			int oldClusterLabel = clusterLabels[i];
			clusterSizes[clusterLabels[i]] --;
			double mintmp = lambda;// cost of starting a new cluster
			int newClusterLabel = nextClusterLabel;
			for(j = 0; j < I + 1; j ++) {
				if(clusterSizes[j] > 0) {
					Dist(i, j) = 0;
					for(k = 0; k < K; k ++) {
						for(s = 0; s < S; s ++) {
							if(Theta(k, i) == s) {
								Dist(i, j) += (1 - W(s * K + k, j)) * (1 - W(s * K + k, j));
							} else {
								Dist(i, j) += W(s * K + k, j) * W(s * K + k, j);
							}
						}
					}
			
					// assign the cluster label as well as the outlier label
					if(Dist(i, j) * lambdaw < mintmp) {
						mintmp = Dist(i, j) * lambdaw;
						newClusterLabel = j;
					}
				}
			}

			clusterLabels[i] = newClusterLabel;
			clusterSizes[newClusterLabel] ++;
			if(mintmp >= lambda) {
				// a new cluster is formed
				// if(J != newClusterLabel) {printf("Error: new state is not J = %d.", J);exit(1);}
				for(s = 0; s < S; s ++) {
					for(k = 0; k < K; k ++) {
						if(Theta(k, i) != s) {
							W(s * K + k, newClusterLabel) = _LOW;
						} else {
							W(s * K + k, newClusterLabel) = 1 - (S - 1) * _LOW;
						}
					}
				}
				for(i1 = 0; i1 < I; i1 ++) {
					Dist(i1, newClusterLabel) = 0;
					for(k = 0; k < K; k ++) {
						if(Theta(k, i) != Theta(k, i1)) {
							Dist(i1, newClusterLabel) += 2;
						}
					}
				}
				J ++;
			}
			if(clusterSizes[oldClusterLabel] == 0) {
				// an old cluster should be removed
				if(false) {
					W(_, oldClusterLabel) = W(_, J - 1);
					Dist(_, oldClusterLabel) = Dist(_, J - 1);
					for(k = 0; k < K * S; k ++) {
						W(k, J - 1) = 0;
						Dist(k, J - 1) = 0;
					}
					for(i1 = 0; i1 < I; i1 ++) {
						if(clusterLabels[i1] == J - 1)
							clusterLabels[i1] = oldClusterLabel;
					}
					clusterSizes[oldClusterLabel] = clusterSizes[J - 1];
					clusterSizes[J - 1] = 0;
				}
				J --;
				nextClusterLabel = oldClusterLabel;
			} else {
				// find the next cluster label
				while(clusterSizes[nextClusterLabel] > 0) {
					nextClusterLabel ++;
					if(nextClusterLabel > I) {
						nextClusterLabel -= I + 1;
					}
				}
			}
			//too many clusters
			if(J >= sqrt(I) * 3) {
				break;
			}
		}

		//loss = ComputeLoss(D, Theta, Y, Mu, Gamma, W, clusterLabels, clusterSizes, lambda, lambdaw);
		//printf("update cluster, Loss function = %3.3f, number of clusters = %d\n", loss, J);
		// test if each cluster can be dismissed
		double dismissLoss = 0, bestLoss = 0;
		int newLabels[I];
		bool canDismiss = false;
		for(j = 0; j < I + 1; j ++) {
			if(clusterSizes[j] > 0) {
				dismissLoss = 0;
				canDismiss = true;
				for(i = 0; i < I; i ++) {
					if(clusterLabels[i] == j) {
						// assign to the next closest cluster
						newLabels[i] = 0;
						bestLoss = 2 * K;
						for(j1 = 0; j1 < I + 1; j1 ++) {
							if(clusterSizes[j1] > 0 && j1 != j && bestLoss > Dist(i, j1)) {
								newLabels[i] = j1;
								bestLoss = Dist(i, j1);
							}
						}
						dismissLoss += bestLoss - Dist(i, j);
						if(dismissLoss * lambdaw > lambda) {
							canDismiss = false;
						}
					}
					if(! canDismiss) {
						break;
					}
				}
				if(canDismiss) {
					// printf("Dismiss cluster %d, %3.3f > %3.3f.\n", j, lambda / lambdaw, dismissLoss);
					for(i = 0; i < I; i ++) {
						if(clusterLabels[i] == j) {
							clusterLabels[i] = newLabels[i];
							clusterSizes[newLabels[i]] ++;
						}
					}
					if(false) {
					if(j < J - 1) {
						// an old cluster should be removed
						W(_, j) = W(_, J - 1);
						Dist(_, j) = Dist(_, J - 1);
						for(i1 = 0; i1 < I; i1 ++) {
							if(clusterLabels[i1] == J - 1)
								clusterLabels[i1] = j;
						}
						clusterSizes[j] = clusterSizes[J - 1];
						j --;// need to recheck this new cluster since it is relabelled
					}
					}
					clusterSizes[j] = 0;
					if(false) {
						for(k = 0; k < K * S; k ++) {
							W(k, J - 1) = 0;
							Dist(k, J - 1) = 0;
						}
						clusterSizes[J - 1] = 0;
					}
					J --;
					nextClusterLabel = j;
				}
			}
		}

		//loss = ComputeLoss(D, Theta, Y, Mu, Gamma, W, clusterLabels, clusterSizes, lambda, lambdaw);
		//printf("dismissed clusters, Loss function = %3.3f, number of clusters = %d\n", loss, J);

		// update mu
		for(n = 0; n < N; n ++) {
			double denom = 0, numer = 0;
			//Initialize Mu as the pooled mean to avoid 0/0
			for(i = 0; i < I; i ++) {
				for(s = 0; s < S; s ++) {
					denom += Gamma(n + s * N, i);
					numer += Y(n, i);
				}
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
			// Initialize Sigma as the pooled variance to avoid 0/0
			denom = 0;
			numer = 0;
			for(i = 0; i < I; i ++) {
				for(s = 0; s < S; s ++) {
					if(Theta(D[n], i) == s) {
						numer += (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s)) * (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s));
						denom ++;
					}
				}
			}
			for(s = 0; s < S; s ++) {
				Sigma(n, s) = numer / denom;
			}
			for(s = 0; s < S; s ++) {
				denom = 0;
				numer = 0;
				for(i = 0; i < I; i ++) {
					if(Theta(D[n], i) == s) {
						numer += (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s)) * (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s));
						denom ++;
					}
				}
				if(denom > 0) {
					Sigma(n, s) = numer / denom;
				}
			}
		}

		// calculate the loss function

		loss = ComputeLoss(D, Theta, Y, Mu, Gamma, W, clusterLabels, clusterSizes, lambda, lambdaw);
		//printf("update distribution, Loss function = %3.3f, number of clusters = %d\n", loss, J);

		if(oldloss - loss < tol * oldloss && itr > 1) {
			break;
		}

		if(J >= sqrt(I) * 3) {
			printf("Warning: too many clusters. Consider using larger lambda.\n");
			break;
		}

		lastLabel = I;
		while(clusterSizes[lastLabel] == 0) {
			lastLabel --;
		}

		if(lastLabel > J * 2) {
			//printf("relabel clusters\n");
			firstLabel = 0, lastLabel = I;
			while(firstLabel < lastLabel) {
				while(clusterSizes[firstLabel] > 0 && firstLabel < I) {
					firstLabel ++;
				}
				while(clusterSizes[lastLabel] == 0 && lastLabel > 0) {
					lastLabel --;
				}
				if(firstLabel < lastLabel) {
					for(i = 0; i < I; i ++) {
						if(clusterLabels[i] == lastLabel) {
							clusterLabels[i] = firstLabel;
						}
					}
					W(_, firstLabel) = W(_, lastLabel);
					Dist(_, firstLabel) = Dist(_, lastLabel);
					clusterSizes[firstLabel] = clusterSizes[lastLabel];
					clusterSizes[lastLabel] = 0;
					nextClusterLabel = lastLabel;
				}
			}
		}
	}

	//printf("relabel clusters\n");
	// Reorder the clusterLabels and the columns of W
	firstLabel = 0, lastLabel = I;
	while(firstLabel < lastLabel) {
		while(clusterSizes[firstLabel] > 0 && firstLabel < I) {
			firstLabel ++;
		}
		while(clusterSizes[lastLabel] == 0 && lastLabel > 0) {
			lastLabel --;
		}
		if(firstLabel < lastLabel) {
			for(i = 0; i < I; i ++) {
				if(clusterLabels[i] == lastLabel) {
					clusterLabels[i] = firstLabel;
				}
			}
			W(_, firstLabel) = W(_, lastLabel);
			for(k = 0; k < K * S; k ++) {
				W(k, lastLabel) = 0;
			}
			clusterSizes[firstLabel] = clusterSizes[lastLabel];
			clusterSizes[lastLabel] = 0;
		}
	}

	Rcpp::List ret = Rcpp::List::create(
					    Rcpp::Named("Theta") = Theta,
					    Rcpp::Named("clusterLabels") = clusterLabels,
					    Rcpp::Named("W") = W,
					    Rcpp::Named("Mu") = Mu,
					    Rcpp::Named("Sigma") = Sigma,
					    Rcpp::Named("loss") = loss,
					    Rcpp::Named("Iter") = itr
					    );
	return( ret );
	
}

double ComputeLoss(IntegerVector D, IntegerMatrix Theta, NumericMatrix Y, NumericMatrix Mu, NumericMatrix Gamma, NumericMatrix W, IntegerVector clusterLabels, IntegerVector clusterSizes, double lambda, double lambdaw) {
	int I = Theta.ncol();
	int K = Theta.nrow();
	int N = D.size();
	int S = Mu.ncol();

	int i, n, k, s, J;
	double loss = 0, tmp;
	for(i = 0; i < I; i ++) {
		//if(clusterLabels[i] >= J) {			J = clusterLabels[i] + 1;		}
		for(n = 0; n < N; n ++) {
			k = D[n];
			s = Theta(k, i);
			loss += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) *
				(Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
		}
	}
	//	printf("loss after data %3.3f\n", loss);
		
	for(i = 0; i < I; i ++) {
		tmp = 0;
		for(k = 0; k < K; k ++) {
			for(s = 0; s < S; s ++) {
				if(Theta(k, i) == s) {
					tmp += (1 - W(s * K + k, clusterLabels[i])) * (1 - W(s * K + k, clusterLabels[i])); 
				} else {
					tmp += W(s * K + k, clusterLabels[i]) * W(s * K + k, clusterLabels[i]);
				}
			}
		}
		//		printf("i = %d, loss = %3.3f, b = %d\n", i, tmp, b[i]);
		loss += lambdaw * tmp;
	}
	//	printf("loss after w %3.3f\n", loss);
	J = 0;
	for(i = 0; i < I + 1; i ++) {
		if(clusterSizes[i] > 0) {
			J ++;
		}
	}
	loss += lambda * (J - 1);
	//	printf("loss for all terms  %3.3f\n", loss);
	return(loss);
}

double ComputeLoss_theta(IntegerVector D, IntegerMatrix Theta, NumericMatrix Y, NumericMatrix Mu, NumericMatrix Gamma) {
	int I = Theta.ncol();
	int N = D.size();

	int i, n, k, s;
	double loss = 0;
	for(i = 0; i < I; i ++) {
		for(n = 0; n < N; n ++) {
			k = D[n];
			s = Theta(k, i);
			loss += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) *
				(Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
		}
	}
	
	return(loss);
}


SEXP RSolveW(SEXP _counts, SEXP _scalar, SEXP _p, SEXP _q) {
	NumericVector counts(_counts);
	double scalar = as<double>(_scalar);
	double p = as<double>(_p);
	double q = as<double>(_q);
	NumericVector ret = SolveW(counts, scalar, p, q);
	return(wrap(ret));
}

//solve for cluster profile
NumericVector SolveW(NumericVector counts, double scalar, double p, double q) {
	int S = counts.size() / 2;
	double _LOW = 1e-6;
	double I1 = 0, I2 = 0;
	for(int s = 0; s < S; s ++) {
		counts[s] ++; // add a pseudo count
		counts[s + S] ++; // add a pseudo count
		I1 += counts[s];
		I2 += counts[s + S];
	}

	//initialize the solvant
	NumericVector sol(S), newsol(S);
	for(int s = 0; s < S; s ++) {
		newsol[s] = (counts[s] + counts[s + S]) / (I1 + I2);
	}

	double d1, d2, exp1, exp2, exp3, exp4;

	int i;
	for(i = 0; i < 100; i ++) {
		//update each variable
		// gradient descent
		
		for(int s = 0; s < S; s ++) {
			sol[s] = newsol[s];
			exp1 = exp(p * log(1 - sol[s]));
			exp2 = exp(p * log(sol[s]));
			exp3 = scalar * exp(q * log(1 - sol[s]));
			exp4 = scalar * exp(q * log(sol[s]));
			
			d1 = - counts[s] * p * exp1 / (1 - sol[s]) + (I1 - counts[s]) * p * exp2 / sol[s]
				- counts[s + S] * q * exp3 / (1 - sol[s]) + (I2 - counts[s + S]) * q * exp4 / sol[s];
			d2 = (counts[s] * exp1 / (1 - sol[s]) / (1 - sol[s]) + (I1 - counts[s]) * exp2 / sol[s] / sol[s]) * p * (p - 1) 
				+ (counts[s + S] * exp3 / (1 - sol[s]) / (1 - sol[s]) + (I2 - counts[s + S]) * exp4 / sol[s] / sol[s]) * q * (q - 1) ;
			newsol[s] -= d1 / d2 / 2;
		}
		
		// projection to the hyperplane
		d1 = 0;
		for(int s = 0; s < S; s ++) {
			d1 += newsol[s];
		}
		for(int s = 0; s < S; s ++) {
			newsol[s] /= d1;
			if(newsol[s] < _LOW) {
				newsol[s] = _LOW;
			} else if(newsol[s] > 1 - _LOW) {
				newsol[s] = 1 - _LOW;
			}
		}

		// check convergence
		d1 = 0;
		for(int s = 0; s < S; s ++) {
			d2 = (newsol[s] - sol[s]) * (newsol[s] - sol[s]);
			if(d1 < d2) {
				d1 = d2;
			}
		}
		if(d1 < 0.0001)
			break;
	}

	//	printf("Number of iterations: %d\n", i);
	
	return(sol);
}


SEXP madbayes_theta(SEXP _Theta, SEXP _Mu, SEXP _D, SEXP _Gamma, SEXP _Y, SEXP _maxitr, SEXP _tol) {
	
	// The following values are 1updated in MCMC iterations
	IntegerMatrix Theta(_Theta); // K by I
	NumericMatrix Mu(_Mu); // N by S
	IntegerVector D(_D); // Length N, valued in {0, 1, ..., K-1}

	// The following is the external information.
	NumericMatrix Gamma(_Gamma); // N*S by I
	NumericMatrix Y(_Y); // N by I

	int maxitr = as<int>(_maxitr);
	double tol = as<double>(_tol);

	// extract the dimensions
	int I = Theta.ncol();
	int S = Mu.ncol();
	int K = Theta.nrow();
	int N = D.size();

	// The following will be computed
	NumericMatrix Sigma(N, S);

	// iterators
	int i, k = 0, s = 0, n, itr;//, likid;
	double loss = 0, oldloss;

	for(itr = 0; itr < maxitr; itr ++) {
		oldloss = loss;

		// update Theta
		for(i = 0; i < I; i ++) {
			for(k = 0; k < K; k ++) {
				double tmp[S];
				for(s = 0; s < S; s ++) {
					// initialize
					tmp[s] = 0;
					for(n = 0; n < N; n ++) {
						if(D[n] == k) {
							tmp[s] += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) * (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
						}
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
		
		//loss = ComputeLoss_theta(D, Theta, Y, Mu, Gamma);
		//	printf("update states, Loss function = %3.3f\n", loss);
		
		// update mu
		for(n = 0; n < N; n ++) {
			double denom = 0, numer = 0;
			//Initialize Mu as the pooled mean to avoid 0/0
			for(i = 0; i < I; i ++) {
				for(s = 0; s < S; s ++) {
					denom += Gamma(n + s * N, i);
					numer += Y(n, i);
				}
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
			// Initialize Sigma as the pooled variance to avoid 0/0
			denom = 0;
			numer = 0;
			for(i = 0; i < I; i ++) {
				for(s = 0; s < S; s ++) {
					if(Theta(D[n], i) == s) {
						numer += (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s)) * (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s));
						denom ++;
					}
				}
			}
			for(s = 0; s < S; s ++) {
				Sigma(n, s) = numer / denom;
			}
			for(s = 0; s < S; s ++) {
				denom = 0;
				numer = 0;
				for(i = 0; i < I; i ++) {
					if(Theta(D[n], i) == s) {
						numer += (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s)) * (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s));
						denom ++;
					}
				}
				if(denom > 0) {
					Sigma(n, s) = numer / denom;
				}
			}
		}
		
		// calculate the loss function
		
		loss = ComputeLoss_theta(D, Theta, Y, Mu, Gamma);
		//printf("initialize distribution, Loss function = %3.3f\n", loss);

		if(itr > 1 && oldloss - loss < tol * oldloss) {
			break;
		}
	}
	
	Rcpp::List ret = Rcpp::List::create(
					    Rcpp::Named("Theta") = Theta,
					    Rcpp::Named("Mu") = Mu,
					    Rcpp::Named("Sigma") = Sigma,
					    Rcpp::Named("loss") = loss,
					    Rcpp::Named("Iter") = itr
					    );
	return( ret );

}

SEXP madbayes_init(SEXP _Theta, SEXP _lambda, SEXP _S, SEXP _maxJ) {
	IntegerMatrix Theta(_Theta);
	double lambda = as<double>(_lambda);
	int S = as<int>(_S);
	int maxJ = as<int>(_maxJ);

	int I = Theta.ncol();
	int K = Theta.nrow();

	//Rcpp::List allClusterLabels(maxJ);
	//NumericVector allLosses(maxJ);

	NumericVector probs(I), newProbs(I), clusterSizes(I + 1), newClusterSizes(I + 1);
	IntegerVector clusterLabels(I), newClusterLabels(I);
	NumericMatrix W(K * S, I + 1), newW(K * S, I + 1);
	double loss, newLoss;
	int centroid;

	// initialize the sampling probabilities
	for(int i = 0; i < I + 1; i ++) {
		for(int k = 0; k < K * S; k ++) {
			W(k, i) = 0;
			newW(k, i) = 0;
		}
		clusterSizes[i] = 0;
		newClusterSizes[i] = 0;
	}

	for(int i = 0; i < I; i ++) {
		for(int k = 0; k < K; k ++) {
			W(Theta(k, i) * K + k, 0) ++;
		}
	}
	for(int k = 0; k < K * S; k ++) {
		W(k, 0) /= (double) I;
	}

	clusterSizes[0] = I;

	loss = 0;
	for(int i = 0; i < I; i ++) {
		probs[i] = 0;
		for(int k = 0; k < K; k ++) {
			for(int s = 0; s < S; s ++) {
				if(Theta(k, i) != s) {
					probs[i] += W(k + K * s, 0) * W(k + K * s, 0);
				} else {
					probs[i] += (1 - W(k + K * s, 0)) * (1 - W(k + K * s, 0));
				}
			}
		}
		clusterLabels[i] = 0;
		loss += probs[i];
	}

	srand(time(NULL));
	
	// Kmeans ++ algorithm to initialize W and cluterLabels
	int J = 1;
	while(J < maxJ) {
		centroid = SampleCentroid(probs);

		// update the clusterLabels for each unit
		// notice that we only need to test if the new center should be the new centroid
		for(int i = 0; i < I; i ++) {
			// compute the distance to the new centroid
			double dist = 0;
			if(i != centroid) {
				for(int k = 0; k < K; k ++) {
					if(Theta(k, i) != Theta(k, centroid)) {
						dist += 2;
					}
				}
			}
			// update the cluster label
			if(probs[i] > dist) {
				newClusterLabels[i] = J;
			}
		}

		// update the W matrix
		for(int j = 0; j <= J; j ++) {
			for(int k = 0; k < K * S; k ++) {
				newW(k, j) = 0;
			}
			newClusterSizes[j] = 0;
			for(int i = 0; i < I; i ++) {
				if(newClusterLabels[i] == j) {
					newClusterSizes[j] ++;
					for(int k = 0; k < K; k ++) {
						newW(k + K * Theta(k, i), j) ++;
					}
				}
			}
			for(int k = 0; k < K; k ++) {
				if(newClusterSizes[j] > 0) {
					for(int s = 0; s < S; s ++) {
						newW(k + K * s, j) /= (double) newClusterSizes[j];
					}
				}
			}
		}

		// remove empty clusters
		int newJ = 0;
		while(true) {
			while(newClusterSizes[J] == 0 && J > 0) {
				J--;
			}
			while(newClusterSizes[newJ] > 0 && newJ < I) {
				newJ ++;
			}
			if(newJ > J) {
				break;
			}
			newW(_, newJ) = newW(_, J);
			newClusterSizes[newJ] = newClusterSizes[J];
			newClusterSizes[J] = 0;
			for(int i = 0; i < I; i ++) {
				if(newClusterLabels[i] == J) {
					newClusterLabels[i] = newJ;
				}
			}
			J--;
		}
		if(J != newJ - 1) {
			printf("Error: unexpected new cluster numbers.\n");
			exit(1);
		}
		
		// compute the new loss
		newLoss = 0;
		for(int i = 0; i < I; i ++) {
			newProbs[i] = 0;
			for(int k = 0; k < K; k ++) {
				for(int s = 0; s < S; s ++) {
					if(Theta(k, i) == s) {
						newProbs[i] += (1 - newW(k + s * K, newClusterLabels[i])) * (1 - newW(k + s * K, newClusterLabels[i]));
					} else {
						newProbs[i] += newW(k + s * K, newClusterLabels[i]) * newW(k + s * K, newClusterLabels[i]);
					}
				}
			}
			newLoss += newProbs[i];
		}
		newLoss += lambda * (J - 1);

		//	printf("new loss = %3.3f\n", newLoss);

		if(newLoss >= loss) {
			break;
		}

		// update
		clusterLabels = newClusterLabels;
		W = newW;
		clusterSizes = newClusterSizes;
		probs = newProbs;
		loss = newLoss;

		//allClusterLabels[J] = Rcpp::clone(clusterLabels);
		//allLosses[J] = loss;

		J = newJ;
		//	printf("update, %d clusters\n", J); 
	}

	Rcpp::List ret = Rcpp::List::create(//Rcpp::Named("allClusterLabels") = allClusterLabels,
					    //Rcpp::Named("allLosses") = allLosses,
					    Rcpp::Named("loss") = loss,
					    Rcpp::Named("clusterLabels") = clusterLabels);
	return(ret);
}

SEXP madbayes_init_kmeanspp(SEXP _Theta, SEXP _S, SEXP _J) {
	IntegerMatrix Theta(_Theta);
	int S = as<int>(_S);
	int J = as<int>(_J);

	int I = Theta.ncol();
	int K = Theta.nrow();

	//Rcpp::List allClusterLabels(J);
	//NumericVector allLosses(J);

	NumericVector probs(I);
	IntegerVector clusterLabels(I);
	IntegerVector centroids(J);

	// initialize the sampling probabilities
	double loss = 0;
	for(int i = 0; i < I; i ++) {
		probs[i] = K * S;
		clusterLabels[i] = 0;
		loss += probs[i];
	}

	srand(time(NULL));
	
	// Kmeans ++ algorithm to initialize W and cluterLabels
	for(int j = 0; j < J; j ++) {
		centroids[j] = SampleCentroid(probs);

		// update the clusterLabels for each unit
		// notice that we only need to test if the new center should be the new centroid
		for(int i = 0; i < I; i ++) {
			// compute the distance to the new centroid
			double dist = 0;
			if(i != centroids[j]) {
				for(int k = 0; k < K; k ++) {
					if(Theta(k, i) != Theta(k, centroids[j])) {
						dist += 2;
					}
				}
			}
			// update the cluster label
			if(probs[i] > dist) {
				loss += dist - probs[i];
				probs[i] = dist;
				clusterLabels[i] = j;
			}
		}
		
		//allClusterLabels[j] = Rcpp::clone(clusterLabels);
		//allLosses[j] = loss;

	}

	Rcpp::List ret = Rcpp::List::create(//Rcpp::Named("allClusterLabels") = allClusterLabels,
					    //Rcpp::Named("allLosses") = allLosses,
					    Rcpp::Named("loss") = loss,
					    Rcpp::Named("clusterLabels") = clusterLabels,
					    Rcpp::Named("centroids") = centroids);
	return(ret);
}

int SampleCentroid(NumericVector probs) {
	int I = probs.size();
	double total = 0;
	int i;
	for(i = 0; i < I; i ++) {
		total += probs[i];
	}
	double rv = ((double) rand() / (RAND_MAX)) * total;
	i = 0;
	total = probs[i];
	while(total < rv) {
		i++;
		total += probs[i];
	}
	return(i);
}
