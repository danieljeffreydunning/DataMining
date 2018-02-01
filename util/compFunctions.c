#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "compFunctions.h"

/*calculate distance to a single point in a 1D array (such as data or query)*/
double pnt2pntDistance(int dim, int query_idx, double *query, int data_idx, double *data) {
	int i;
	double distance = 0.0;

	for (i = 0; i < dim; i++) {
		distance += pow(data[data_idx+i] - query[query_idx+i], 2);
	}
    distance = sqrt(distance);

	return distance;
}

/*calculate distance to a centroid in a 2D array (cluster_centroid) */
double pnt2centDistance(int dim, int clust_idx, int query_idx, double *query, double **cluster_centroid) {
	int i;
	double distance = 0.0;

	for (i = 0; i < dim; i++) {
		distance += pow(cluster_centroid[clust_idx][i] - query[query_idx+i], 2);
	}
    distance = sqrt(distance);

	return distance;
}

/*calculate the dot product between a data point and a 2D array of values. returns int because we need the floor*/
int dotprod(int dim, double *data, double **r, int d_idx, int r_idx) {
	int i, i_sum;
    double sum = 0;

	for (i = 0; i < dim; i++) {
		sum += data[d_idx+i] * r[r_idx][i];
	}
    i_sum = (int) sum;

	return i_sum;
}

/*calculate distance to kdtree cluster boundry. Returns 0 if it is within the boundries*/
double pnt2bdry(int dim, double **cluster_bdry, double *query, int qidx, int clust_idx) {
    int i;
    double q0, min_bdry, max_bdry, calcdist = 0;

    for (i = 0; i < 2 * dim; i+=2) { //for each dimension boundry
        q0 = query[qidx+(i/2)];
        min_bdry = cluster_bdry[clust_idx][i];
        max_bdry = cluster_bdry[clust_idx][i+1];

        if (q0 < min_bdry) //smaller than the min boundry
            calcdist += pow(q0 - min_bdry, 2);
        else if (q0 > max_bdry) //larger than the max boundry
            calcdist += pow(q0 - max_bdry, 2);
        else //right in the middle, nothing changes
            calcdist += 0;
    }

    calcdist = sqrt(calcdist);

    return calcdist;
}

