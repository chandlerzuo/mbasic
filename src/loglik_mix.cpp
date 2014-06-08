#include "loglik_mix.h"

SEXP loglik_mix( SEXP _W, SEXP _p, SEXP _zeta, SEXP _probz, SEXP _PDF, SEXP _V ) {
	NumericMatrix PDF( _PDF );
        NumericMatrix W( _W );
	NumericVector p( _p );
	NumericMatrix V( _V );
	double zeta = as<double>( _zeta );
	NumericVector probz( _probz );

	// extract the dimensions
	int I = p.size();
	int K = PDF.nrow();
	int J = probz.size();
	int S = PDF.ncol() / I;

	double loglik = 0;
	int i, j, k, s;

	double dbl_min = exp( - 64 );

	for( i = 0; i < I; i ++ ){
		double exp1 = log( zeta );
		for( k = 0; k < K; k ++ ){
			double maxexp = PDF( k, i );
			for( s = 0; s < S; s ++ ){
				if( maxexp < PDF( k, i + I * s ) )
					maxexp = PDF( k, i + I * s );
			}
			double tmp = 0; 
			for( s = 0; s < S; s ++ ){
				if( s == S - 1 )
					tmp += p( i ) * exp( PDF( k, i ) - maxexp );
				else
					tmp += ( 1 - p( i ) ) * V( k, s ) * exp( PDF( k, i + (s+1) * I ) - maxexp );					
			}
			if( tmp < dbl_min )
				tmp = dbl_min;
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
				for( s = 0; s < S; s ++ ){
					if( s == S - 1 )
						tmp += W( k, j ) * exp( PDF( k, i ) - maxexp );
					else
						tmp += ( 1 - W( k, j ) ) * V( k, s ) * exp( PDF( k, i + ( s + 1 ) * I ) - maxexp );
				}
				if( tmp < dbl_min )
					tmp = dbl_min;
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


