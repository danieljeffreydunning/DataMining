#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 16
#define K 3 //at least 2 
#define DIM 2
#define Q 10

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define LOG2(X) log((X)) / log(2)

//function prototypes
int kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, int *cluster_assign);
void calculateCentroids(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign);
void assignData(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign);
int search_kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double *query, double *result_pt);
double calculateClusterDistance(int dim, int data_idx, double *data, int clust_idx, double **cluster_centroid);
void initializeCentroids(int dim, int ndata, double *data, int k, double **cluster_centroid, int start_idx, int m);
//end prototypes


void calculateCentroids(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign) {
	//keep in mind cluster_start values are in terms of cluster_assign
	int i, j, l, loopstart, loopend;
	double temp_centroid[DIM] = {0.0}; 
	
	for (i = 0; i < k; i++) {
		loopstart = cluster_start[i]*dim;
		loopend =  loopstart + cluster_size[i]*dim;
		//printf("loop start %d loop end %d \n", loopstart, loopend);
		//for each centroid, getting the sum of all points, each of their dimensions, in temp_centroid
		for (j = loopstart; j < loopend; j+=dim) {
			for (l = 0; l < dim; l++) {	
				temp_centroid[l] += data[j+l];
			}
		} 

		//now get the averages for each dimension in temp_centroid and put it in cluster_centroid
		//printf("cluster centroid ");
		for (j = 0; j < dim; j++) {
			//printf("%f\t,", temp_centroid[j]); 
			cluster_centroid[i][j] = temp_centroid[j] / k;
			temp_centroid[j] = 0.0;
		}
		printf("\n");
	}	


}


double calculateClusterDistance(int dim, int data_idx, double *data, int clust_idx, double **cluster_centroid) {
	int i;
	double distance = 0.0;

	for (i = 0; i < dim; i++) {
		distance += pow(cluster_centroid[clust_idx][i] - data[data_idx+i], 2);
	}

	distance = sqrt(distance);

	return distance;
}

void assignData(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign) {
	int i, j, l, cent_idx, temp_data_idx = 0, clust_size_idx_cnt = 0;
	double temp_distance, distance = 10000.0;
	double *temp_data;

	temp_data = (double *)malloc(sizeof(double) * ndata * dim);

	//intitialize temp data
	for ( i = 0; i < ndata*dim; i++) {
		temp_data[i] = data[i];
	}


	for (i = 0; i < ndata*dim; i+=dim) { //each data point
		for (j = 0; j < k; j++) { //each centroid 
			temp_distance = calculateClusterDistance(dim, i, data, j, cluster_centroid);

			distance = MIN(distance, temp_distance);

			if (distance == temp_distance)
				cent_idx = j;
		}

		cluster_assign[i/dim] = cent_idx;
	}

	/*printf("\n");
	for (i = 0; i < ndata; i++) {
		printf("\t%d,", cluster_assign[i]);
	}
	printf("\n");*/

	for (l = 0; l < k; l++) { // for each cluster 
		for (i = 0; i < ndata*dim; i+=dim) { //for each data point	
			//rearrange data set
			if (cluster_assign[i/dim] == l) {
				cluster_size[l]++; //keep track of each cluster size
				for (j = 0; j < dim; j++) {
					data[temp_data_idx+j] = temp_data[i+j];
				}	
			}
		}
	}

	//get cluster_starts from sizes
	for (l = 0; l < k; l++) {
		cluster_start[l] = clust_size_idx_cnt;
		clust_size_idx_cnt += cluster_size[l];
	}

	free(temp_data);
}

void initializeCentroids(int dim, int ndata, double *data, int k, double **cluster_centroid, int start_idx, int m) {
	int i, j, l, centroid_point;
	double distance = 0.0, point_min = 10000.0, total_max = 0.0;
	//start_idx will be 0 in initialization case, but will be the cluster to be redone if calling this function for an empty cluster
	

	for (j = 0; j < ndata*dim; j+=dim) { //for each data point
		for (i = start_idx; i < m; i++) { //for given m centroids
			for (l = 0; l < dim; l++) {
				distance += pow(cluster_centroid[i][l] - data[j+l], 2); 		
			}	
			distance = sqrt(distance);
			//get min distance to a centroid for that point
			point_min = MIN(point_min, distance);
			//printf("min distance for point %d is %f\n", j/dim, point_min);
			 
		}
		//get max min distance for each data point
		total_max = MAX(total_max, point_min);
	
		//if we own the current max min, this will currently be the centroid point
		if (total_max == point_min) {
			centroid_point = j;	
		}

		//reset point_min
		point_min = 10000.0;

	}
		
	//now put the new centroid in its place
	for (i = 0; i < dim; i++) {
		cluster_centroid[m][i] = data[centroid_point+i];
	}
}


int main() {
	
	int i, j, kcheck = 1, rint;
	double r;
	double *data, **cluster_centroid, *query, *result;
	int *cluster_assign, *cluster_size,  *cluster_start;

	data = (double *)malloc(sizeof(double) * N * DIM);
	cluster_assign = (int *)malloc(sizeof(int) * N);
	cluster_size = (int *)malloc(sizeof(int) * K);
	cluster_start = (int *)malloc(sizeof(int) * K);
	cluster_centroid = (double **)malloc(sizeof(double *) * K);

	query = (double *)malloc(sizeof(double) * Q * DIM);
	result = (double *)malloc(sizeof(double) * DIM);

	//initialize 2d arrays
	for (i = 0; i < K; i++) {
		cluster_centroid[i] = (double *)malloc(sizeof(double) * DIM);
	}

	//printf("cluster assign at i: ");
	//initialize cluster_assign
	for (i = 0; i < N; i++) {
		cluster_assign[i] = i * DIM;
	//	 printf("%d," , cluster_assign[i]);
	}
	
	//printf("\n");


	srand(time(NULL));
	
	for (i = 0; i < N*DIM; i+=DIM) {
		for (j = 0; j < DIM; j++) {
			r = ((double)rand() / (double)(RAND_MAX)) * 100.0;
			data[i+j] = r;
			//printf("%f,", r);
		}
		//printf("\n");
	}
	printf("\n");

	//initialize centroids
	rint = rand() % N;
	for (i = 0; i < DIM; i++) {
		cluster_centroid[0][i] = data[rint+i];
	}

	for (i = 0; i < K; i++) {
		cluster_size[i] = 0;
	}
/*	for (i = 0; i < N*DIM; i+=DIM) {
		td = calculateDistance(DIM, i, data, cluster_centroid); 
		d = MAX(d, td);
		if (td == d) {
			cent_idx = i;
		}
	}
	for (i = 0; i < dim; i++) {
		
	}*/
	while (kcheck < K) {
		initializeCentroids(DIM, N, data, K, cluster_centroid, 0, kcheck);
		kcheck++;
	}

	

	/*for (i = 0; i < Q * DIM; i++) {
		r = ((double)rand() / (double)(RAND_MAX)) * 100.0;
		query[i] = r; 
		//printf("%f,", r);
	}*/


	//printf("dank\n");
	printf("\n");
	for (i = 0; i < N*DIM; i+=DIM) {
		printf("%d) ", i/DIM);
		for (j = 0; j < DIM; j++) {
			printf("%f,", data[i+j]);
		}
		printf("\n");
	}

	printf("\nCentroid\n");
	for (i = 0; i < K; i++) {
		for (j = 0; j < DIM; j++) {
			printf("%f,", cluster_centroid[i][j]);
		}
		printf("\n");
	}


	assignData(DIM, N, data, K, cluster_size, cluster_start, cluster_centroid, cluster_assign);
	calculateCentroids(DIM, N, data, K, cluster_size, cluster_start, cluster_centroid, cluster_assign);

	printf("\nCentroid\n");
	for (i = 0; i < K; i++) {
		for (j = 0; j < DIM; j++) {
			printf("%f,", cluster_centroid[i][j]);
		}
		printf("\n");
	}
	/*printf("\nQuery: ");
	for (i = 0; i < Q*DIM; i++) {
		printf("%f,", query[i]);	
	} 
	printf("\n\n");*/


	/*avg_checks = search_kdtree(DIM, N, data, K, cluster_size, cluster_start, cluster_boundry, query, result);

	printf("Average number of checks was %d\n", avg_checks);
	printf("\n");*/

	return 0;
}

