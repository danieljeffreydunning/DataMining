#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 16
#define K 2 //even power of 2 
#define DIM 2
#define Q 10

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define LOG2(X) log((X)) / log(2)





int main() {
	
	int i, j;
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
		cluster_centroid[i] = (double *)malloc(sizeof(double) * DIM * 2);
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

