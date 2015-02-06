#include "e_step.h"

SEXP e_step( SEXP _W, SEXP _P, SEXP _zeta, SEXP _probz, SEXP _PDF ) {

	// the following parameters are inputs and are not updated
	NumericMatrix PDF( _PDF );

	// the following parameters serve both as input, but will be updated in M step as output
        NumericMatrix W( _W );
	NumericMatrix P( _P );
	double zeta = as<double>( _zeta );
	NumericVector probz( _probz );
	
	// extract the dimensions
	int I = P.nrow();
	int S = P.ncol();
	int K = W.nrow();
	int J = probz.size();

	double _LOW = 1e-10;

	// the following quantities are outputs
	NumericMatrix Theta_mean( K, I * S );
	NumericVector b_mean( I );
	NumericMatrix Theta_b( I, S );
	NumericVector Z_mean( J );
	NumericMatrix W_max( K, J*S );
	NumericMatrix predZ( I, J );

	// iterators
	int i, j, k, s;//, likid;

	// Intermediate matrices
	NumericMatrix FP( K, I );
	NumericMatrix FW( K, I * J );
	NumericMatrix Zcond( I, J );
	
	for( k = 0; k < K; k ++ ){
		for( i = 0; i < I; i ++ ){
			FP( k, i ) = 0; 
			double tmpexp[ S ];
			for( s = 0; s < S; s ++ ){
				if( P( i, s ) < _LOW )
					P( i, s ) = _LOW;
				tmpexp[ s ] = PDF( k, i + s * I ) + log( P( i, s ) );
			}
			double maxexp = tmpexp[ 0 ];
			for( s = 1; s < S; s ++ )
				if( maxexp < tmpexp[ s ] )
					maxexp = tmpexp[ s ];
			for( s = 0; s < S; s ++ )
				FP( k, i ) += exp( tmpexp[ s ] - maxexp );
			FP( k, i ) = log( FP( k, i ) ) +  maxexp;
			
			for( j = 0; j < J; j ++ ){
				FW( k, i + j * I ) = 0;
				for( s = 0; s < S; s ++ ){
					if( W( k, j + J * s ) < _LOW )
						W( k, j + J * s ) = _LOW;
					tmpexp[ s ] = PDF( k, i + I * s ) + log( W( k, j + J * s ) );
				}
				maxexp = tmpexp[ 0 ];
				for( s = 1; s < S; s ++ )
					if( maxexp < tmpexp[ s ] )
						maxexp = tmpexp[ s ];
				for( s = 0; s < S; s ++ )
					FW( k, i + j * I ) += exp( tmpexp[ s ] - maxexp );
				FW( k, i + j * I ) = log( FW( k, i + j * I ) ) + maxexp;
			}
		}
	}

	/*
	// compute old likelihood
	double oldlik = 0;
	for( i = 0; i < I; i ++ ){
		double exp1 = log( zeta );
		for( k = 0; k < K; k ++ ){
			exp1 += FP( k, i );
		}
		double exp2 = 0;
		double tmpexp[ J ];
		for( j = 0; j < J; j ++ ){
			tmpexp[ j ] = log( probz[ j ] );
			for( k = 0; k < K; k ++ )
				tmpexp[ j ] +=  FW( k, i + j * I );
		}
		double maxexp = tmpexp[ 0 ];
		for( j = 1; j < J; j ++ )
			if( maxexp < tmpexp[ j ] )
				maxexp = tmpexp[ j ];
		for( j = 0; j < J; j ++ )
			exp2 += exp( tmpexp[ j ] - maxexp );
		exp2 += maxexp;
		exp2 += log( 1 - zeta );
		if( exp1 > exp2 )
			oldlik += log( 1 + exp( exp2 - exp1 ) ) + exp1;
		else
			oldlik += log( 1 + exp( exp1 - exp2 ) ) + exp2;
	}
	*/
	
	for( j = 0; j < J; j ++ )
		Z_mean( j ) = 0;

	// compute b, Z
	for( i = 0; i < I; i ++ ){

		double exp1,exp0;
		// numerator
		exp1 = log( zeta );
		for( k = 0; k < K; k ++ )
			exp1 += FP( k, i );

		// denominator
		double tmpexp[ J ];
		for( j = 0; j < J; j ++ ){
			tmpexp[ j ] = log( 1-zeta ) + log( probz( j ) );
			for( k = 0; k < K; k ++ )
				tmpexp[ j ] += FW( k, i + I * j );
		}
		double maxexp = tmpexp[0];
		for( j = 1; j < J; j ++ )
			if( tmpexp[ j ] > maxexp )
				maxexp = tmpexp[ j ];
		exp0 = 0;
		for( j = 0; j < J; j ++ )
			exp0 += exp( tmpexp[ j ] - maxexp );
		exp0 = log( exp0 );
		exp0 += maxexp;

		// this matrix will be recycled
		for( j = 0; j < J; j ++ )
			Zcond( i, j ) = exp( tmpexp[ j ] - exp0 );

		//calculate b
		if( exp1 > exp0 )
			b_mean( i ) = 1 / ( 1 + exp( exp0 - exp1 ) );
		else 
			b_mean( i ) = 1 - 1 / ( 1 + exp( exp1 - exp0 ) );

		//calculate Z
		double tmpe[J];
		for( j = 0; j < J; j ++ ){
			double tmp = exp1 + log( probz[ j ] );
			if( tmp > tmpexp[ j ] )
				tmpe[ j ] = tmp + log( 1 + exp( tmpexp[j] - tmp ) );
			else
				tmpe[ j ] = tmpexp[ j ] + log( 1 + exp( tmp - tmpexp[j] ) );
		}
		
		maxexp = tmpe[ 0 ];
		for( j = 1; j < J; j ++ )
			if( maxexp < tmpe[ j ] )
				maxexp = tmpe[ j ];
		double denom = 0;
		for( j = 0; j < J; j ++ ){
			tmpe[ j ] -= maxexp;
			denom += exp( tmpe[ j ] );
		}
		denom = log( denom );

		for( j = 0; j < J; j ++ ){
			predZ( i, j ) = exp( tmpe[ j ] - denom );
			Z_mean( j ) += predZ( i, j );
		}

	}


	for( i = 0; i < I; i ++ )
		for( s = 0; s < S; s ++ )
			Theta_b( i, s ) = 0;

	for( k = 0; k < K; k ++ )
		for( j = 0; j < J; j ++ )
			for( s = 0; s < S; s ++ )
				W_max( k, j + J * s ) = 0;

	for( i = 0; i < I; i ++ ){
		for( k = 0; k < K; k ++ ){
			for( s = 0; s < S; s ++ ){
				double tmp1;
				tmp1 = b_mean( i ) * exp( PDF( k, i + s * I ) - FP( k, i ) ) * P( i, s );
				Theta_b( i, s ) += tmp1;
				if(i == 64) {
					printf("i=%d, k=%d, s=%d, Theta_b(i, s) += %3.3f\n b_mean=%3.3f, PDF=%3.3f, FP=%3.3f, P=%3.3f\n\n", i, k, s, tmp1, b_mean(i), PDF(k, i + s * I), FP(k, i), P(i, s));
				} 
				Theta_mean( k, i + I * s ) = tmp1;
				for( j = 0; j < J; j ++ ){
					double tmp2 = ( 1 - b_mean( i ) ) * Zcond( i, j ) * exp( PDF( k, i + I * s )  - FW( k, i + I * j ) ) *  W( k, j + J * s );
					Theta_mean( k, i + I * s ) += tmp2;
					W_max( k, j + J * s ) += tmp2;
				}
				
				if( Theta_mean( k, i + I * s ) < _LOW )
					Theta_mean( k, i + I * s ) = _LOW;
				else if( Theta_mean( k, i + I * s ) > 1 - _LOW )
					Theta_mean( k, i + I * s ) = 1 - _LOW;
				
			}
		}
		
	}
	
	/*
	for( k = 0; k < K; k ++ ){
		for( j = 0; j < J; j ++ ){
			double total = 0;
			for( s = 0; s < S; s ++ )
				total += W_max( k, j + J * s );
			total = 1;
			for( s = 0; s < S; s ++ ){
				if( W_max( k, j + J * s ) < _LOW * total )
					W_max( k, j + J * s ) = _LOW;
				else if( W_max( k, j + J * s ) > ( 1 -_LOW ) * total )
					W_max( k, j + J * s ) = 1 - _LOW;
				else if( total < _LOW )
					W_max( k, j + J * s ) = 1 / S;
				else
					W_max( k, j + J * s ) /= total;
			}
		}
		}*/

	zeta = 0;
	for( i = 0; i < I; i ++ )
		zeta += b_mean[ i ];
	zeta /= I;
	if( zeta < _LOW )
		zeta = _LOW;
	else if( zeta > 1 - _LOW )
		zeta = 1 - _LOW;
	
	for( j = 0; j < J; j ++ ){
		Z_mean( j ) /= I;
		if( Z_mean( j ) < _LOW )
			Z_mean( j ) = _LOW;
		else if( Z_mean( j ) > 1 - _LOW )
			Z_mean( j ) = 1 - _LOW;
	}


	for( i = 0; i < I; i ++ ){
		double total = 0;
		for( s = 0; s < S; s ++ )
			total += Theta_b( i, s );
		if(i == 64)
			printf("total = %3.20f\t", total);
		for( s = 0; s < S; s ++ ){
			if(i == 64)
				printf("Theta_b(%d, %d) = %3.20\t", i, s, Theta_b(i, s)); 
			if( Theta_b( i, s ) < _LOW * total )
				Theta_b( i, s ) = _LOW;
			else if( Theta_b( i, s ) > ( 1-_LOW ) * total )
				Theta_b( i, s ) = 1 - _LOW;
			else if( total < _LOW )
				Theta_b( i, s ) = 1 / S;
			else
				Theta_b( i, s ) /= total;
			if(i == 64)
				printf("--> %3.20\n", Theta_b(i, s)); 
		}
	}

	Rcpp::List ret = Rcpp::List::create(
					    Rcpp::Named("P") = Theta_b,
					    Rcpp::Named("zeta") = zeta,
					    Rcpp::Named("probz") = Z_mean,
					    Rcpp::Named("W") = W_max,
					    Rcpp::Named("Theta_mean") = Theta_mean,
					    Rcpp::Named("b_prob") = b_mean,
					    Rcpp::Named("predZ") = predZ
					    //,Rcpp::Named("oldlik") = oldlik
					    );
	
	return( ret );
	
}
