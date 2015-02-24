#include "e_step.h"

SEXP e_step(SEXP _W, SEXP _P, SEXP _V, SEXP _zeta, SEXP _probz, SEXP _PDF, SEXP _factor, SEXP _statemap) {

	// the following parameters are inputs and are not updated
	NumericMatrix PDF(_PDF); // NM by I
	IntegerVector factor(_factor); // N by K
	IntegerVector statemap(_statemap); // M by S

	// the following parameters serve both as input, but will be updated in M step as output
        NumericMatrix W(_W);
	NumericMatrix P(_P);
	NumericMatrix V(_V); // N by M
	double zeta = as<double>(_zeta);
	NumericVector probz(_probz);
	
	// extract the dimensions
	int S = P.ncol();
	int I = P.nrow();
	int J = probz.size();
	int K = W.nrow() / S;
	int N = factor.size();
	int M = statemap.size();

	double _LOW = 1e-10;

	// the following quantities are outputs
	NumericMatrix Theta_mean(K * S, I);
	NumericVector b_mean(I);
	NumericMatrix Theta_b(I, S);
	NumericVector Z_mean(J);
	NumericMatrix W_max(K * S, J);
	NumericMatrix V_max(N, M);
	NumericMatrix V_norm(N, S);
	NumericMatrix predZ(I, J);
	NumericMatrix Theta_nu(N * M, I);

	// iterators
	int i, j, k, s, m, n;//, likid;

	// Intermediate matrices
	NumericMatrix jointPDF(K * S, I);
	NumericMatrix FP(K, I);
	NumericMatrix FW(K * J, I);
	NumericMatrix FV(N * S, I);
	NumericMatrix Zcond(I, J);
	
	//compute joint PDF
	for(i = 0; i < I; i ++) {
		for(k = 0; k < K; k ++) {
			for(s = 0; s < S; s ++) {
				jointPDF(k + K * s, i) = 0;
				for(n = 0; n < N; n ++) {
					if(factor[n] == k) {
						int m1 = 0;
						while(statemap[m1] != s) {
							m1 ++;
						}
						double maxexp = PDF(n + N * m1, i);
						for(m = m1; m < M; m ++) {
							if(statemap[m] == s && PDF(n + N * m, i) > maxexp) {
								maxexp = PDF(n + N * m, i);
							}
						}
						double tmp = 0;
						for(m = 0; m < M; m ++) {
							if(statemap[m] == s) {
								tmp += V(n, m) * exp(PDF(n + N * m, i) - maxexp);
							}
						}
						if(tmp < _LOW) {
							tmp = _LOW;
						}
						jointPDF(k + K * s, i) += log(tmp) + maxexp;
					}
				}
			}
		}
	}

	for(k = 0; k < K; k ++){
		for(i = 0; i < I; i ++){
			// compute FP
			double tmpexp[S], maxexp;
			if(zeta > 0) {
				FP(k, i) = 0; 
				for(s = 0; s < S; s ++){
					if(P(i, s) < _LOW)
						P(i, s) = _LOW;
					tmpexp[s] = jointPDF(k + K * s, i) + log(P(i, s));
				}
				maxexp = tmpexp[0];
				for(s = 1; s < S; s ++)
					if(maxexp < tmpexp[s])
						maxexp = tmpexp[s];
				for(s = 0; s < S; s ++)
					FP(k, i) += exp(tmpexp[s] - maxexp);
				FP(k, i) = log(FP(k, i)) +  maxexp;
			}
			
			// compute FW
			for(j = 0; j < J; j ++){
				FW(k + K * j, i) = 0;
				for(s = 0; s < S; s ++){
					if(W(k + K * s, j) < _LOW)
						W(k + K * s, j) = _LOW;
					tmpexp[s] = jointPDF(k + K * s, i) + log(W(k + K * s, j));
				}
				maxexp = tmpexp[0];
				for(s = 1; s < S; s ++)
					if(maxexp < tmpexp[s])
						maxexp = tmpexp[s];
				for(s = 0; s < S; s ++)
					FW(k + K * j, i) += exp(tmpexp[s] - maxexp);
				FW(k + K * j, i) = log(FW(k + K * j, i)) + maxexp;
			}
		}
	}

	//compute FV
	for(i = 0; i < I; i ++) {
		for(n = 0; n < N; n ++) {
			for(s = 0; s < S; s ++) {
				int m1 = 0;
				while(statemap[m1] != s) {
					m1 ++;
				}
				double maxexp = PDF(n + N * m1, i);
				for(m = m1; m < M; m ++) {
					if(statemap[m] == s && PDF(n + N * m, i) > maxexp) {
						maxexp = PDF(n + N * m, i);
					}
				}
				double tmp = 0;
				for(m = 0; m < M; m ++) {
					if(statemap[m] == s) {
						tmp += exp(PDF(n + N * m, i) - maxexp) * V(n, m);
					}
				}
				if(tmp < _LOW) {
					tmp = _LOW;
				}
				FV(n + N * s, i) = log(tmp) + maxexp;
			}
		}
	}
							   
	for(j = 0; j < J; j ++)
		Z_mean(j) = 0;

	// compute b, Z
	for(i = 0; i < I; i ++){

		double exp1 = 0,exp0 = 0;
		// numerator
		if(zeta > 0) {
			exp1 = log(zeta);
			for(k = 0; k < K; k ++)
				exp1 += FP(k, i);
		}

		// denominator
		double tmpexp[J];
		for(j = 0; j < J; j ++){
			tmpexp[j] = log(1-zeta) + log(probz(j));
			for(k = 0; k < K; k ++)
				tmpexp[j] += FW(k + K * j, i);
		}
		double maxexp = tmpexp[0];
		for(j = 1; j < J; j ++)
			if(tmpexp[j] > maxexp)
				maxexp = tmpexp[j];
		exp0 = 0;
		for(j = 0; j < J; j ++)
			exp0 += exp(tmpexp[j] - maxexp);
		exp0 = log(exp0);
		exp0 += maxexp;

		// this matrix will be recycled
		for(j = 0; j < J; j ++)
			Zcond(i, j) = exp(tmpexp[j] - exp0);

		//calculate b
		if(zeta > 0) {
			if(exp1 > exp0)
				b_mean(i) = 1 / (1 + exp(exp0 - exp1));
			else 
				b_mean(i) = 1 - 1 / (1 + exp(exp1 - exp0));
		} else {
			b_mean(i) = 0;
		}

		//calculate Z
		if(zeta > 0) {
			double tmpe[J];
			for(j = 0; j < J; j ++){
				double tmp = exp1 + log(probz[j]);
				if(tmp > tmpexp[j])
					tmpe[j] = tmp + log(1 + exp(tmpexp[j] - tmp));
				else
					tmpe[j] = tmpexp[j] + log(1 + exp(tmp - tmpexp[j]));
			}
			
			maxexp = tmpe[0];
			for(j = 1; j < J; j ++)
				if(maxexp < tmpe[j])
					maxexp = tmpe[j];
			double denom = 0;
			for(j = 0; j < J; j ++){
				tmpe[j] -= maxexp;
				denom += exp(tmpe[j]);
			}
			denom = log(denom);
			
			for(j = 0; j < J; j ++){
				predZ(i, j) = exp(tmpe[j] - denom);
				Z_mean(j) += predZ(i, j);
			}
		} else {
			for(j = 0; j < J; j ++){
				predZ(i, j) = Zcond(i, j);
				Z_mean(j) += predZ(i, j);
			}
		}

	}

	for(i = 0; i < I; i ++)
		for(s = 0; s < S; s ++)
			Theta_b(i, s) = 0;

	for(k = 0; k < K; k ++)
		for(j = 0; j < J; j ++)
			for(s = 0; s < S; s ++)
				W_max(k + K * s, j) = 0;

	for(n = 0; n < N; n ++) {
		for(m = 0; m < M; m ++) {
			V_max(n, m) = 0;
			for(i = 0; i < I; i ++)
				Theta_nu(n + N * m, i) = 0;
		}
		for(s = 0; s < S; s ++) {
			V_norm(n, s) = 0;
		}
	}
	

	for(i = 0; i < I; i ++){
		for(k = 0; k < K; k ++){
			for(s = 0; s < S; s ++){
				double tmp1;
				tmp1 = b_mean(i) * exp(jointPDF(k + K * s, i) - FP(k, i)) * P(i, s);
				Theta_b(i, s) += tmp1;
				Theta_mean(k + s * K, i) = tmp1;
				for(j = 0; j < J; j ++){
					double tmp2 = (1 - b_mean(i)) * Zcond(i, j) * exp(jointPDF(k + K * s, i)  - FW(k + K * j, i)) *  W(k + K * s, j);
					Theta_mean(k + K * s, i) += tmp2;
					W_max(k + K * s, j) += tmp2;
				}
				
				if(Theta_mean(k + K * s, i) < _LOW)
					Theta_mean(k + K * s, i) = _LOW;
				else if(Theta_mean(k + K * s, i) > 1 - _LOW)
					Theta_mean(k + K * s, i) = 1 - _LOW;
				
				for(n = 0; n < N; n ++) {
					if(factor[n] == k) {
						for(m = 0; m < M; m ++) {
							if(statemap[m] == s) {
								Theta_nu(n + N * m, i) = Theta_mean(k + K * s, i) * V(n, m) * exp(PDF(n + N * m, i) - FV(n + N * s, i));
								if(Theta_nu(n + N * m, i) < _LOW)
									Theta_nu(n + N * m, i) = _LOW;
								else if(Theta_nu(n + N * m, i) > 1 - _LOW)
									Theta_nu(n + N * m, i) = 1 - _LOW;
								V_max(n, m) += Theta_nu(n + N * m, i);
								V_norm(n, s) += Theta_nu(n + N * m, i);
							}
						}
					}
				}

			}
		}
	}
	
	if(zeta > 0) {
		zeta = 0;
		for(i = 0; i < I; i ++)
			zeta += b_mean[i];
		zeta /= I;
		if(zeta < _LOW)
			zeta = _LOW;
		else if(zeta > 1 - _LOW)
			zeta = 1 - _LOW;
	}
	
	for(j = 0; j < J; j ++){
		Z_mean(j) /= I;
		if(Z_mean(j) < _LOW)
			Z_mean(j) = _LOW;
		else if(Z_mean(j) > 1 - _LOW)
			Z_mean(j) = 1 - _LOW;
	}


	for(i = 0; i < I; i ++){
		double total = 0;
		for(s = 0; s < S; s ++)
			total += Theta_b(i, s);
		for(s = 0; s < S; s ++){
			if(Theta_b(i, s) < _LOW * total)
				Theta_b(i, s) = _LOW;
			else if(Theta_b(i, s) > (1-_LOW) * total)
				Theta_b(i, s) = 1 - _LOW;
			else if(total < _LOW)
				Theta_b(i, s) = 1 / (double)S;
			else
				Theta_b(i, s) /= total;
		}
	}
	
	for(n = 0; n < N; n ++) {
		for(s = 0; s < S; s ++) {
			for(m = 0; m < M; m ++) {
				if(statemap[m] == s) {
					V_max(n, m) /= V_norm(n, s);
				}
			}
		}
	}

	Rcpp::List ret = Rcpp::List::create(
					    Rcpp::Named("P") = Theta_b,
					    Rcpp::Named("zeta") = zeta,
					    Rcpp::Named("probz") = Z_mean,
					    Rcpp::Named("W") = W_max,
					    Rcpp::Named("Theta") = Theta_mean,
					    Rcpp::Named("Theta_nu") = Theta_nu,
					    Rcpp::Named("V") = V_max,
					    Rcpp::Named("b_prob") = b_mean,
					    Rcpp::Named("predZ") = predZ
					   );
	
	return(ret);
	
}
