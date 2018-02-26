#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "bisecting_kmeans.h"
#include "util/compFunctions.h"


void calculateInitCentroids(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int clust1, int clust2) {
	//keep in mind cluster_start values are in terms of cluster_assign
	int i, j, l, loopstart, loopend, range;
	double *temp_centroid; 

    temp_centroid = (double *)malloc(sizeof(double) * dim);

    for (i = 0; i < dim; i++) {
        temp_centroid[i] = 0.0;
    }
	
	range = clust2 - clust1;

	for (i = clust1; i < clust2+1; i+=range) {//for each cluster
		loopstart = cluster_start[i];
		loopend =  loopstart + cluster_size[i]*dim;
		//printf("loop start %d loop end %d \n", loopstart, loopend);
		//for each centroid, getting the sum of all points, each of their dimensions, in temp_centroid
		for (j = loopstart; j < loopend; j+=dim) {
			for (l = 0; l < dim; l++) {	
				temp_centroid[l] += data[j+l];
			}
		} 

		//now get the averages for each dimension in temp_centroid and put it in cluster_centroid
		//printf("temp cluster centroid ");
		for (j = 0; j < dim; j++) {
			//printf("%f\t,", temp_centroid[j]); 
			cluster_centroid[i][j] = temp_centroid[j] / cluster_size[i];
			temp_centroid[j] = 0.0;
		}
		//printf("\n");
	}	

    free(temp_centroid);
}

int assignInitData(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int clust1, int clust2) { //ndata does not refer to whole data set, only to amount in the specified cluster
	int i, j, l, cent_idx, temp_data_idx = cluster_start[clust1], clust_size_idx_cnt = cluster_start[clust1], empty_idx = k+5, range, temp_size = 0;
	int clust1_size = 0, clust2_size = 0;
	double temp_distance, distance = DBL_MAX;
	double *temp_data, *temp_assign;

	temp_data = (double *)malloc(sizeof(double) * ndata * dim);
	temp_assign = (double *)malloc(sizeof(double) * ndata);

	range = clust2 - clust1;

	//intitialize temp data
	for (i = cluster_start[clust1]; i < cluster_start[clust1]+ndata*dim; i++) {
		temp_data[i] = data[i];
	}

	for (i = cluster_start[clust1]; i < cluster_start[clust1]+cluster_size[clust1]*dim; i+=dim) { //each data point
		for (j = clust1; j < clust2+1; j+=range) { //each centroid 
			//temp_distance = calculateClusterDistance(dim, i, data, j, cluster_centroid);
            temp_distance = pnt2centDistance(dim, j, i, data, cluster_centroid);
			distance = MIN(distance, temp_distance);

			if (distance == temp_distance)
				cent_idx = j;
		}
		distance = DBL_MAX;
		cluster_assign[i/dim] = cent_idx;
	}

	for (l = clust1; l < clust2+1; l+=range) { // for each cluster	
		for (i =cluster_start[clust1]; i < cluster_start[clust1]+ndata*dim; i+=dim) { //for each data point in the cluster original cluster	
			//rearrange data set
			if (cluster_assign[i/dim] == l) {
				temp_size++; //keep track of each cluster size
				for (j = 0; j < dim; j++) {
					data[temp_data_idx+j] = temp_data[i+j];
				}	
				temp_assign[temp_data_idx/dim] = l;
				temp_data_idx+=dim;
			}
		}
		cluster_size[l] = temp_size;
		temp_size = 0;
	}

	//rearrange cluster_assign
	//printf("\nCluster Assign ");
	for (i = cluster_start[clust1]/dim; i < cluster_start[clust1]/dim+ndata; i++) {
		cluster_assign[i] = temp_assign[i];
		//printf("%d,", cluster_assign[i]);
	}
	//printf("\n");


	//get cluster_starts from sizes
	//printf("cluster size, ");
	for (l = clust1; l < clust2+1; l+=range) {
		cluster_start[l] = clust_size_idx_cnt;
		clust_size_idx_cnt += cluster_size[l]*dim;
		//printf("%d,", cluster_size[l]); 
	}
	//printf("\n");

	for (i = 0; i < k; i++) {
		if (cluster_size[i] == 0) {
			empty_idx = i;
			break;
		}	
	}

	free(temp_assign);
	free(temp_data);

	return empty_idx;
}

int bkInitCentroids(int dim, int ndata, double *data, int *cluster_size, int *cluster_start, double **cluster_centroid, int place_idx, int cent_idx, int cent_count, int *cluster_assign) {
	int i, j, l, centroid_point, empty_idx, sums_idx;
	double distance = 0.0, total_max = 0.0, sums_hold = 0.0, sums_max = 0.0;
	double *sums_arr;
	//start_idx will be 0 in initialization case, but will be the cluster to be redone if calling this function for an empty cluster
	sums_arr = (double *)malloc(sizeof(double) * cent_count);

	for (i = 0; i < cent_count; i++) {
		sums_arr[i] = 0.0;
	}

	for (j = cluster_start[cent_idx]; j < cluster_start[cent_idx]+cluster_size[cent_idx]*dim; j+=dim) { //for each data point in the cluster
		for (l = 0; l < dim; l++) {
				distance += pow(cluster_centroid[cent_idx][l] - data[j+l], 2); 		
			}	
			distance = sqrt(distance);
			//printf("min distance for point %d is %f\n", j/dim, point_min);

		//get max distance for each data point
		total_max = MAX(total_max, distance);
	
		//if we own the current max min, this will currently be the centroid point
		if (total_max == distance) {
			centroid_point = j;	
		}
		distance = 0.0;
	}
		
	//now put the new centroid in its place
	for (i = 0; i < dim; i++) {
		cluster_centroid[place_idx][i] = data[centroid_point+i];
	}


	empty_idx = assignInitData(dim, cluster_size[cent_idx], data, 2, cluster_size, cluster_start, cluster_centroid, cluster_assign, cent_idx, place_idx);
	calculateInitCentroids(dim, cluster_size[cent_idx], data, 2, cluster_size, cluster_start, cluster_centroid, cluster_assign, cent_idx, place_idx);

	for (i = 0; i < cent_count; i++) { //for each cluster
		for (l = cluster_start[i]; l < cluster_start[i]+cluster_size[i]*dim; l+=dim) {
			for (j = 0; j < dim; j++) {
				sums_hold += pow(cluster_centroid[i][j] - data[l+j], 2); 
			}
			sums_hold = sums_hold; // for large data sets, don't want to overflow our array
			sums_arr[i] += sums_hold;
			sums_hold = 0.0;
		}
	}

	for (i = 0; i < cent_count; i++) {
		sums_hold = sums_arr[i];
		sums_max = MAX(sums_max, sums_hold);
		if (sums_max == sums_hold) {
			sums_idx = i;
		}
	}

	/*printf("\nCentroid round %d\n", place_idx);
	for (i = 0; i < cent_count; i++) {
		for (j = 0; j < dim; j++) {
			printf("%f,", cluster_centroid[i][j]);
		}
		printf("\n");
	}

	printf("\nsums idx %d\n", sums_idx);*/

	free(sums_arr);

	return sums_idx;
}

void runBKmeans(int dim, int k, int ndata, double *data, double **cluster_centroid) {
	
	int i, j, kcheck = 1, rint, cycles, pointcnt, next_clust = 0;
	double r;
	int *cluster_assign, *cluster_size,  *cluster_start;

	cluster_assign = (int *)malloc(sizeof(int) * ndata);
	cluster_size = (int *)malloc(sizeof(int) * k);
	cluster_start = (int *)malloc(sizeof(int) * k);


	//printf("cluster assign at i: ");
	//initialize cluster_assign
	for (i = 0; i < ndata; i++) {
		cluster_assign[i] = 0;
	//	 printf("%d," , cluster_assign[i]);
	}
	
	//printf("\n");


	srand(23);
	

	//initialize centroids
	rint = rand() % ndata;
	for (i = 0; i < dim; i++) {
		cluster_centroid[0][i] = data[rint+i];
	}

	for (i = 0; i < k; i++) {
		cluster_size[i] = 0;
		cluster_start[i] = ndata*dim;
	}
	
	cluster_size[0] = ndata;
	cluster_start[0] = 0;

/*	for (i = 0; i < N*DIM; i+=DIM) {
		td = calculateDistance(DIM, i, data, cluster_centroid); 
		d = MAX(d, td);
		if (td == d) {
			cent_idx = i;
		}
	}
	for (i = 0; i < dim; i++) {
		
	}*/

	/*printf("\n");
	for (i = 0; i < N*DIM; i+=DIM) {
		printf("%d) ", i/DIM);
		for (j = 0; j < DIM; j++) {
			printf("%f,", data[i+j]);
		}
		printf("\n");
	}*/

	/*printf("\nCentroid\n");
	for (i = 0; i < 1; i++) {
		for (j = 0; j < DIM; j++) {
			printf("%f,", cluster_centroid[i][j]);
		}
		printf("\n");
	}*/

	while (kcheck < k) {
		next_clust = bkInitCentroids(dim, ndata, data, cluster_size, cluster_start, cluster_centroid, kcheck, next_clust, kcheck+1, cluster_assign);
		kcheck++;		
	}

	/*for (i = 0; i < Q * DIM; i++) {
		r = ((double)rand() / (double)(RAND_MAX)) * 100.0;
		query[i] = r; 
		printf("%f,", r);
	}*/


	//printf("dank\n");

	//cycles = kmeans(DIM, N, data, K, cluster_size, cluster_start, cluster_radius, cluster_centroid, cluster_assign);

	/*printf("\nNew Centroid\n");
	for (i = 0; i < K; i++) {
		for (j = 0; j < DIM; j++) {
			printf("%f,", cluster_centroid[i][j]);
		}
		printf("\n");
	}*/

	/*printf("\nCentroid sizes ");
	for (i = 0; i < K; i++) {
		printf("%d,", cluster_size[i]);
	}
	printf("\n");*/

	//printf("\ncycles %d\n", cycles);

	//printf("\n");

	/*pointcnt = search_kmeans(DIM, N, data, K, cluster_size, cluster_start, cluster_radius, cluster_centroid, query, result);
	for (i = 0; i < Q; i++) {
		for (j = 0; j < DIM; j++) {
			printf("%f,", result[i*DIM+j]);
		}
		printf("\n");
	}
	printf("number of points checked %d\n", pointcnt/Q);
*/
	/*printf("\nQuery: ");
	for (i = 0; i < Q*DIM; i++) {
		printf("%f,", query[i]);	
	} 
	printf("\n\n");*/


	/*avg_checks = search_kdtree(DIM, N, data, K, cluster_size, cluster_start, cluster_boundry, query, result);
	printf("Average number of checks was %d\n", avg_checks);
	printf("\n");*/

	free(cluster_start);
	free(cluster_size);
}

