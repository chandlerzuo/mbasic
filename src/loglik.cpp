#include "loglik.h"

SEXP loglik( SEXP _W, SEXP _P, SEXP _zeta, SEXP _probz, SEXP _PDF ) {
	NumericMatrix PDF( _PDF );
        NumericMatrix W( _W );
	NumericMatrix P( _P );
	double zeta = as<double>( _zeta );
	NumericVector probz( _probz );

	// extract the dimensions
	int I = P.nrow();
	int K = PDF.nrow();
	int J = probz.size();
	int S = P.ncol();

	double loglik = 0;
	int i, j, k, s;

	for( i = 0; i < I; i ++ ){
		double exp1 = log( zeta );
		for( k = 0; k < K; k ++ ){
			double maxexp = PDF( k, i );
			for( s = 0; s < S; s ++ ){
				if( maxexp < PDF( k, i + I * s ) )
					maxexp = PDF( k, i + I * s );
			}
			double tmp = 0; 
			for( s = 0; s < S; s ++ )
				tmp += P( i, s ) * exp( PDF( k, i + s * I ) - maxexp );
			exp1 += maxexp + log( tmp );
		}

		double tmp_exp[J];
		for( j = 0; j < J; j ++ ){
			tmp_exp[j] = log( probz[j] );// bug found to cause decreasing likelihood 
			for( k = 0; k < K; k ++ ){
				double maxexp = PDF( k, i );
				for( s = 0; s < S; s ++ ){
					if( maxexp < PDF( k, i + I * s ) )
						maxexp = PDF( k, i + I * s );
				}
				double tmp = 0; 
				for( s = 0; s < S; s ++ )
					tmp += W( k, j + J * s ) * exp( PDF( k, i + s * I ) - maxexp );
				tmp_exp[ j ] += maxexp + log( tmp );
			}
		}
		double max_exp = tmp_exp[0];
		for( j = 1; j < J; j ++ )
			if( tmp_exp[j] > max_exp )
				max_exp = tmp_exp[j];

		double exp2 = 0;
		for( j = 0; j < J; j ++ )
			exp2 += exp( tmp_exp[j] - max_exp );

		exp2 = log( exp2 ) + max_exp;
		exp2 += log( 1 - zeta );
		
		if( exp1 > exp2 )
			loglik += exp1 + log( 1 + exp( exp2 - exp1 ) );
		else
			loglik += exp2 + log( 1 + exp( exp1 - exp2 ) );
	}

	return( wrap( loglik ) );
}


