#include "loglik_theta.h"

SEXP loglik_theta(SEXP _W, SEXP _P, SEXP _zeta, SEXP _probz, SEXP _Theta) {
	NumericMatrix Theta(_Theta);
        NumericMatrix W(_W);
	NumericMatrix P(_P);
	double zeta = as<double>(_zeta);
	NumericVector probz(_probz);

	// extract the dimensions
	int I = P.nrow();
	int J = probz.size();
	int S = P.ncol();
	int K = W.nrow() / S;

	double loglik = 0;
	int i, j, k, s;

	
	for(i = 0; i < I; i ++){
		double exp1 = log(zeta);
		double tmpexp[ J ];
		
		for(j = 0; j < J; j ++)
			tmpexp[ j ] = log(probz[ j ]);
		for(k = 0; k < K; k ++){
			for(s = 0; s < S; s ++){
				if(Theta(k, i + I * s) > 0){
					exp1 += log(P(i, s));
					for(j = 0; j < J; j ++)
						tmpexp[ j ] += log(W(k + K * s, j));
				}
			}
		}

		double maxexp = tmpexp[ 0 ];
		for(j = 1; j < J; j ++)
			if(maxexp < tmpexp[ j ])
				maxexp = tmpexp[ j ];

		double exp2 = 0;
		for(j = 0; j < J; j ++)
			exp2 += exp(tmpexp[ j ] - maxexp);
		exp2 = log(exp2) + maxexp + log(1 - zeta);

		if(exp1 > exp2)
			loglik += exp1 + log(1 + exp(exp2 - exp1));
		else
			loglik += exp2 + log(1 + exp(exp1 - exp2));

	}

	return(wrap(loglik));
}


