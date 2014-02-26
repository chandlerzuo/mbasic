#include "avg_score.h"

SEXP avg_score( SEXP _gene, SEXP _map, SEXP _gc, SEXP _binmap, SEXP _bingc ) {
	NumericMatrix gene( _gene );
        NumericVector map( _map );
	NumericVector gc( _gc );
	double binmap = as<double>( _binmap );
	double bingc = as<double>( _bingc );

	int i,j,k;
	int n = gene.nrow();
	NumericMatrix avemgc(n,2);

	for( i = 0; i < n; i ++ ){
		k=0;
		for( j = (int)((gene(i,0)+1)/binmap); j <= (int)((gene(i,1)+1)/binmap); j ++ ){
			k++;
			avemgc(i,0) += map[j];
		}
		avemgc(i,0) /= k;
	}

	for( i = 0; i < n; i ++ ){
		k=0;
		for( j = (int)((gene(i,0)+1)/bingc); j <= (int)((gene(i,1)+1)/bingc); j ++ ){
			k++;
			avemgc(i,1) += gc[j];
		}
		avemgc(i,1) /= k;
	}

	return( wrap( avemgc ) );
}


