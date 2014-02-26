#include "hamming.h"

SEXP hamming( SEXP _X ) {
	NumericMatrix X( _X );

	// extract the dimensions
	int I = X.ncol();
	int K = X.nrow();

	NumericMatrix D( I, I );

	for( int i = 0; i < I; i ++ ){
		for( int j = i + 1; j < I; j ++ ){
			for( int k = 0; k < K; k ++ )
				if( X( k, i ) != X( k, j ) )
					D( i, j ) ++;
			
		}
	}

	return( wrap( D ) );
}


