#include "silhouette.h"

// compute the silhouette score

SEXP silhouette( SEXP _Theta, SEXP _States, SEXP _J) {

	IntegerMatrix Theta(_Theta);
	IntegerVector States(_States);
	int J = as<int>(_J);

	int I = States.size();
	int K = Theta.nrow();

	int i, j, i1, k;

	double mean_sil = 0;

	NumericMatrix Dist(I, I);
	
	// compute the distance matrix
	for(i = 0; i < I; i ++) {
		Dist(i, i) = 0;
		for(i1 = i + 1; i1 < I; i1 ++) {
			Dist(i, i1) = 0;
			for(k = 0; k < K; k ++) {
				if(Theta(k, i) != Theta(k, i1)) {
					Dist(i, i1) ++;
				}
			}
			Dist(i1, i) = Dist(i, i1);
		}
	}

	// compute the size of each cluster
	double n_cluster[J];
	for(j = 0; j < J; j ++)
		n_cluster[j] = 0;
	for(i = 0; i < I; i ++) {
		n_cluster[States[i]] ++;
	}

	for(i = 0; i < I; i ++) {
		// compute the silhouette score for one point
		// compute the distance to each cluster
		double mean_dist[J];
		for(j = 0; j < J; j ++)
			mean_dist[j] = 0;
		for(i1 = 0; i1 < I; i1 ++) {
			mean_dist[States[i1]] += Dist(i, i1);
		}

		// compute a_i: the average distance to other point in the same cluster
		
		if(n_cluster[States[i]] <= 1) {
			// if the cluster has only one element
			mean_sil --;
			continue;
		}

		double a = mean_dist[States[i]] / (n_cluster[States[i]] - 1);
		
		//compute b_i: the minimum distance to the closest cluster that i does not belong to

		//step 1: find the closest cluster
		for(j = 0; j < J; j ++) {
			if(j != States[i])
				mean_dist[j] /= n_cluster[j];
		}

		double b;
		int minJ = 0;
		if(States[i] == 0)
			minJ = 1;
		
		b = mean_dist[minJ];

		for(j = 0; j < J; j ++) {
			if(mean_dist[j] < b && j != States[i]) {
				minJ = j;
				b = mean_dist[j];
			}
		}

		// add to overall silhouette score

		if(a > b)
			mean_sil += b / a - 1;
		else if(a < b)
			mean_sil += 1 - a / b;
		//		printf("%d\t%3.3f\t%3.3f\t%3.3f\t%d\n", i, a, b, mean_sil, minJ);

	}

	return(wrap(mean_sil / I));
}
