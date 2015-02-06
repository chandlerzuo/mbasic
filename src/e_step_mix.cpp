#include "e_step_mix.h"

SEXP e_step_mix(SEXP _W, SEXP _p, SEXP _zeta, SEXP _probz, SEXP _PDF, SEXP _V) {

	// the following parameters are inputs and are not updated
	NumericMatrix PDF(_PDF);

	// the following parameters serve both as input, but will be updated in M step as output
        NumericMatrix W(_W);
	NumericVector p(_p);
	NumericMatrix V(_V);
	double zeta = as<double>(_zeta);
	NumericVector probz(_probz);
	
	// extract the dimensions
	int I = p.size();
	int S = PDF.ncol() / I;
	int K = W.nrow();
	int J = probz.size();
	double _LOW = 1e-10;

	// the following quantities are outputs
	NumericMatrix Theta_mean(K, I * S);
	NumericVector b_mean(I);
	NumericVector Theta_b(I);
	NumericVector Theta_b_total(I);
	NumericVector Z_mean(J);
	NumericMatrix W_max(K, J);
	NumericMatrix W_total(K, J);
	NumericMatrix V_max(K, S - 1);
	NumericMatrix predZ(I, J);

	// iterators
	int i, j, k, s;//, likid;

	// Intermediate matrices
	NumericMatrix FP(K, I);
	NumericMatrix FW(K, I * J);
	NumericMatrix Zcond(I, J);
	
	for(k = 0; k < K; k ++){
		for(i = 0; i < I; i ++){
			FP(k, i) = 0; 
			double tmpexp[ S ];
			if(p(i) < _LOW)
				p(i) = _LOW;
			tmpexp[ 0 ] = PDF(k, i) + log(p(i));
			for(s = 0; s < S - 1; s ++){
				if((1 - p(i) ) * V(k, s) < _LOW)
					tmpexp[ s+1 ] = PDF(k, i + (s + 1) * I) + log(_LOW);
				else
					tmpexp[ s+1 ] = PDF(k, i + (s + 1) * I) + log((1 - p(i)) * V(k, s));
			}
			double maxexp = tmpexp[ 0 ];
			for(s = 1; s < S; s ++)
				if(maxexp < tmpexp[ s ])
					maxexp = tmpexp[ s ];
			for(s = 0; s < S; s ++)
				FP(k, i) += exp(tmpexp[ s ] - maxexp);
			FP(k, i) = log(FP(k, i)) +  maxexp;
			
			for(j = 0; j < J; j ++){
				FW(k, i + j * I) = 0;
				if(W(k, j) < _LOW)
					W(k, j) = _LOW;
				tmpexp[ 0 ] = PDF(k, i) + log(W(k, j));
				for(s = 0; s < S - 1; s ++){
					if((1 - W(k, j)) * V(k, s) < _LOW)
						tmpexp[ s+1 ] = PDF(k, i + I * (s + 1)) + log(_LOW);
					else
						tmpexp[ s+1 ] = PDF(k, i + I * (1 + s)) + log((1 - W(k, j)) * V(k, s));							}
				maxexp = tmpexp[ 0 ];
				for(s = 1; s < S; s ++)
					if(maxexp < tmpexp[ s ])
						maxexp = tmpexp[ s ];
				for(s = 0; s < S; s ++)
					FW(k, i + j * I) += exp(tmpexp[ s ] - maxexp);
				FW(k, i + j * I) = log(FW(k, i + j * I)) + maxexp;
			}
		}
	}

	for(j = 0; j < J; j ++)
		Z_mean(j) = 0;

	// compute b, Z
	for(i = 0; i < I; i ++){

		double exp1,exp0;
		// numerator
		exp1 = log(zeta);
		for(k = 0; k < K; k ++)
			exp1 += FP(k, i);

		// denominator
		double tmpexp[ J ];
		for(j = 0; j < J; j ++){
			tmpexp[ j ] = log(1-zeta) + log(probz(j));
			for(k = 0; k < K; k ++)
				tmpexp[ j ] += FW(k, i + I * j);
		}
		double maxexp = tmpexp[0];
		for(j = 1; j < J; j ++)
			if(tmpexp[ j ] > maxexp)
				maxexp = tmpexp[ j ];
		exp0 = 0;
		for(j = 0; j < J; j ++)
			exp0 += exp(tmpexp[ j ] - maxexp);
		exp0 = log(exp0);
		exp0 += maxexp;

		// this matrix will be recycled
		for(j = 0; j < J; j ++)
			Zcond(i, j) = exp(tmpexp[ j ] - exp0);

		//calculate b
		if(exp1 > exp0)
			b_mean(i) = 1 / (1 + exp(exp0 - exp1));
		else 
			b_mean(i) = 1 - 1 / (1 + exp(exp1 - exp0));

		//calculate Z
		double tmpe[J];
		for(j = 0; j < J; j ++){
			double tmp = exp1 + log(probz[ j ]);
			if(tmp > tmpexp[ j ])
				tmpe[ j ] = tmp + log(1 + exp(tmpexp[j] - tmp));
			else
				tmpe[ j ] = tmpexp[ j ] + log(1 + exp(tmp - tmpexp[j]));
		}
		
		maxexp = tmpe[ 0 ];
		for(j = 1; j < J; j ++)
			if(maxexp < tmpe[ j ])
				maxexp = tmpe[ j ];
		double denom = 0;
		for(j = 0; j < J; j ++){
			tmpe[ j ] -= maxexp;
			denom += exp(tmpe[ j ]);
		}
		denom = log(denom);

		for(j = 0; j < J; j ++){
			predZ(i, j) = exp(tmpe[ j ] - denom);
			Z_mean(j) += predZ(i, j);
		}

	}

	for(i = 0; i < I; i ++){
		Theta_b(i) = 0;
		Theta_b_total(i) = 0;
	}

	for(k = 0; k < K; k ++)
		for(j = 0; j < J; j ++){
			W_max(k, j) = 0;
			W_total(k, j) = 0;
		}

	for(i = 0; i < I; i ++){
		for(k = 0; k < K; k ++){
			for(s = 0; s < S; s ++){
				double tmp1;
				if(s == S-1){
					tmp1 = b_mean(i) * exp(PDF(k, i) - FP(k, i)) * p(i);
					Theta_b(i) += tmp1;
					Theta_mean(k, i) = tmp1;
				} else {
					tmp1 = b_mean(i) * exp(PDF(k, i + I * (s + 1)) - FP(k, i)) * (1 - p(i)) * V(k, s);
					V_max(k, s) += tmp1;
					Theta_mean(k, i + I * (s + 1)) = tmp1;
				}
				Theta_b_total(i) += tmp1;
				for(j = 0; j < J; j ++){
					double tmp2;
					if(s == S - 1){
						tmp2 = (1 - b_mean(i)) * Zcond(i, j) * exp(PDF(k, i)  - FW(k, i + I * j)) *  W(k, j);
						W_max(k, j) += tmp2;
						Theta_mean(k, i) += tmp2;
					} else {
						tmp2 = (1 - b_mean(i)) * Zcond(i, j) * exp(PDF(k, i + I * (s + 1))  - FW(k, i + I * j)) *  (1 - W(k, j)) * V(k, s);
						V_max(k, s) += tmp2;
						Theta_mean(k, i + I * (s + 1)) += tmp2;
					}
					W_total(k, j) += tmp2;
				}
				if(Theta_mean(k, i + I * s) < _LOW)
					Theta_mean(k, i + I * s) = _LOW;
				else if(Theta_mean(k, i + I * s) > 1 - _LOW)
					Theta_mean(k, i + I * s) = 1 - _LOW;
				
			}
		}
	}

	for(i = 0; i < I; i ++)
		zeta += b_mean[ i ];
	zeta /= I;
	if(zeta < _LOW)
		zeta = _LOW;
	else if(zeta > 1 - _LOW)
		zeta = 1 - _LOW;
	
	for(j = 0; j < J; j ++){
		Z_mean(j) /= I;
		if(Z_mean(j) < _LOW)
			Z_mean(j) = _LOW;
		else if(Z_mean(j) > 1 - _LOW)
			Z_mean(j) = 1 - _LOW;
	}

	for(k = 0; k < K; k ++)
		for(j = 0; j < J; j ++)
			if(W_max(k, j) < _LOW * W_total(k, j))
				W_max(k, j) = _LOW;
			else if(W_max(k, j) > (1 - _LOW) * W_total(k, j))
				W_max(k, j) = 1 - _LOW;
			else
				W_max(k, j) /= W_total(k, j);


	for(k = 0; k < K; k ++){
		double total = 0;
		for(s = 0; s < S - 1; s ++)
			total += V_max(k, s);
		for(s = 0; s < S - 1; s ++)
			if(V_max(k, s) < _LOW * total)
				V_max(k, s) = _LOW;
			else if(V_max(k, s) > (1 - _LOW) * total)
				V_max(k, s) = 1 - _LOW;
			else if(total < _LOW)
				V_max(k, s) = 1 / (double)(S - 1);
			else
				V_max(k, s) /= total;
	}

	for(i = 0; i < I; i ++){
		if(Theta_b(i) < _LOW * Theta_b_total(i))
			Theta_b(i) = _LOW;
		else if(Theta_b(i) > (1-_LOW) * Theta_b_total(i))
			Theta_b(i) = 1 - _LOW;
		else if(Theta_b_total(i) < _LOW)
			Theta_b(i) = 1 / 2;
		else
			Theta_b(i) /= Theta_b_total(i);
	}

	Rcpp::List ret = Rcpp::List::create(
					    Rcpp::Named("p") = Theta_b,
					    Rcpp::Named("zeta") = zeta,
					    Rcpp::Named("probz") = Z_mean,
					    Rcpp::Named("W") = W_max,
					    Rcpp::Named("Theta_mean") = Theta_mean,
					    Rcpp::Named("b_prob") = b_mean,
					    Rcpp::Named("predZ") = predZ,
					    Rcpp::Named("V") = V_max
					    //,Rcpp::Named("oldlik") = oldlik
					   );
	
	return(ret);
	
}
