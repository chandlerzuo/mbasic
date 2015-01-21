#include "e_step_theta.h"

SEXP e_step_theta1(SEXP _W, SEXP _P, SEXP _probz, SEXP _Theta) {

	// the following parameters are inputs and are not updated
	NumericMatrix Theta(_Theta);

	// the following parameters serve both as input, but will be updated in M step as output
        NumericMatrix W(_W);
	NumericMatrix P(_P);
	NumericVector probz(_probz);
	
	// extract the dimensions
	int I = P.nrow();
	int S = P.ncol();
	int K = W.nrow();
	int J = probz.size();

	double _LOW = 1e-10;

	// the following quantities are outputs
	NumericVector b_mean(I);
	NumericVector Z_mean(J);
	NumericMatrix W_max(K, J*S);
	NumericMatrix predZ(I, J);

	// iterators
	int i, j, k, s;//, likid;

	// Intermediate matrices
	NumericVector TP(I);
	NumericMatrix TW(I, J);
	NumericMatrix Zcond(I, J);
	
	for(i = 0; i < I; i ++){
		TP(i) = 0;
		for(j = 0; j < J; j ++)
			TW(i, j) = 0;
		for(k = 0; k < K; k ++){
			for(s = 0; s < S; s ++){
				if(Theta(k, i + I * s) > 0){
					TP(i) += log(P(i, s));
					for(j = 0; j < J; j ++)
						TW(i, j) += log(W(k, j + J * s));
				}
			}
		}
	}

	for(k = 0; k < K; k ++)
		for(j = 0; j < J; j ++)
			for(s = 0; s < S; s ++)
				W_max(k, j + s * J) = 0;
	
	for(j = 0; j < J; j ++)
		Z_mean(j) = 0;

	for(i = 0; i < I; i ++){
		// b_mean
		b_mean(i) = 0;

		// predZ

		double tmpexp[ J ];
		for(j = 0; j < J; j ++){
			tmpexp[ j ] = log(probz[ j ]) + TW(i, j);
		}

		double maxexp = tmpexp[ 0 ];
		for(j = 1; j < J; j ++)
			if(maxexp < tmpexp[ j ])
				maxexp = tmpexp[ j ];

		double total = 0;
		for(j = 0; j < J; j ++)
			total += exp(tmpexp[ j ] - maxexp);

		for(j = 0; j < J; j ++){
			predZ(i, j) = exp(tmpexp[ j ] - maxexp) / total;
			Z_mean(j) += predZ(i, j);
		}

		// Zcond
		for(j = 0; j < J; j ++){
			double exp1 = predZ(i, j);
			//exp2 = 1 - b_mean(i);
			if(exp1 < _LOW)
				Zcond(i, j) = _LOW;
			else if(exp1 >= (1 - _LOW))
				Zcond(i, j) = 1 - _LOW;
			else
				Zcond(i, j) = exp1;

			for(k = 0; k < K; k ++)
				for(s = 0; s < S; s ++){
					if(Theta(k, i + I * s) > 0){
						W_max(k, j + J * s) += Zcond(i, j);
						//	printf("%d, %d\t%d\t%d\t%f\t%f\t%f\t%f\n", i, k, j, s, W_max(k, j + J * s), Zcond(i, j), predZ(i, j), b_mean(i), probz(j));
						break;
					}
				}
		}

	}

	for(j = 0; j < J; j ++){
		Z_mean(j) /= I;
		if(Z_mean(j) < _LOW)
			Z_mean(j) = _LOW;
		else if(Z_mean(j) > 1 - _LOW)
			Z_mean(j) = 1 - _LOW;
	}

	for(k = 0; k < K; k ++)
		for(j = 0; j < J; j ++){
			double total = 0;
			for(s = 0; s < S; s ++)
				total += W_max(k, j + J * s);
			for(s = 0; s < S; s ++){
				if(total == 0)
					W_max(k, j + s * J) = 1 / S;
				else if(W_max(k, j + s * J) < _LOW * total)
					W_max(k, j + s * J) = _LOW;
				else if(W_max(k, j + s * J) > (1 - _LOW) * total)
					W_max(k, j + s * J) = 1 - _LOW;
				else 
					W_max(k, j + s * J) /= total;
			}
		}
	
	Rcpp::List ret = Rcpp::List::create(
					    Rcpp::Named("zeta") = 0,
					    Rcpp::Named("probz") = Z_mean,
					    Rcpp::Named("W") = W_max,
					    Rcpp::Named("b_prob") = b_mean,
					    Rcpp::Named("predZ") = predZ
					    //,Rcpp::Named("oldlik") = oldlik
					   );
	
	return(ret);
	
}
