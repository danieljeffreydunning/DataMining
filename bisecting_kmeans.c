#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 64
#define K 8 //at least 2 
#define DIM 2
#define Q 1
#define SIZEFACTOR 10

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define LOG2(X) log((X)) / log(2)

//function prototypes
int kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, int *cluster_assign);
void calculateCentroids(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign);
void calculateInitCentroids(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int clust1, int clust2);
int assignData(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign);
int assignInitData(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int clust1, int clust2);
int search_kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, double *query, double *result_pt);
double calculateClusterDistance(int dim, int data_idx, double *data, int clust_idx, double **cluster_centroid);
double calculatePointDistance(int dim, int clust_idx, int query_idx, double *query, double **cluster_centroid);
double distance(int dim, int query_idx, double *query, int data_idx, double *data);
int initializeCentroids(int dim, int ndata, double *data, int *cluster_size, int *cluster_start, double **cluster_centroid, int place_idx, int cent_idx, int cent_count, int *cluster_assign);
//end prototypes

int search_kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, double *query, double *result_pt) {
	int i, j, l, min_clust_idx, min_point_idx, count = 0;
	double cent_dist = 10000.0, point_dist = 10000.0, temp_dist, rad_comp_dist;
	double *cent_comp_arr;

	cent_comp_arr = (double *)malloc(sizeof(double) * k);

	for (j = 0; j < Q; j++) { // for each query
		//initialize comp arr
		for (i = 0; i < k; i++) {
			cent_comp_arr[i] = 0.0;
		}
		for (i = 0; i < k; i++) { //for each cluster
			temp_dist = calculatePointDistance(dim, i, j, query, cluster_centroid); //calculate distance to each centroid

			cent_comp_arr[i] = temp_dist;
			
			cent_dist = MIN(cent_dist, temp_dist);
			if (cent_dist == temp_dist) {
				min_clust_idx = i;
			}

		}
		printf("closest cluster %d\n", min_clust_idx);
		
		//calculate distance to each point in the closest cluster
		for (l = cluster_start[min_clust_idx]; l < cluster_start[min_clust_idx]+cluster_size[min_clust_idx]*dim; l+=dim) {
			temp_dist = distance(dim, j, query, l, data);
			count++;

			point_dist = MIN(point_dist, temp_dist);
			if (point_dist == temp_dist) {
				min_point_idx = l;
			}	
		}
		//check if any cluster radius is closer than the min distance. If so, check points in that cluster
		for (i = 0; i < k; i++) {
			if (i != min_clust_idx) { //maker sure its not the closest cluster we already did
				rad_comp_dist = cent_comp_arr[i] - cluster_radius[i]; //distance to the edge of that cluster
				if (rad_comp_dist <= point_dist) {
					//calculate distance to each point in the next closest cluster(s)
					for (l = cluster_start[i]; l < cluster_start[i]+cluster_size[i]*dim; l+=dim) {
						temp_dist = distance(dim, j, query, l, data);
						count++;

						point_dist = MIN(point_dist, temp_dist);
						if (point_dist == temp_dist) {
							min_point_idx = l;
						}	
					}
	
				}
			}
		}

		for (i = 0; i < dim; i++) {
			result_pt[j*dim+i] = data[min_point_idx+i];
		}

	}

	return count;
}


/*int kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, int *cluster_assign) {
	int i, j, doneflag = 0, cycle_cnt = 1, loopstart, loopend, stallflag = 0, empty_idx, cycle_check = 0, f = 0;
	int *temp_assign, *cycle_check_assign, max_dist = 0.0, temp_max;
	
	temp_assign = (int *)malloc(sizeof(int) * ndata);
	cycle_check_assign = (int *)malloc(sizeof(int) * ndata);

	//first call of assignData
	empty_idx = assignData(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign);

	//make sure there are now empty clusters
	while ((empty_idx < k) && (f < k)) {
		f++;
		initializeCentroids(dim, ndata, data, k, cluster_centroid, empty_idx, k, k); 
		empty_idx = assignData(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign);
	}	
	f = 0;

	calculateCentroids(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign);

	while ((!doneflag) && (stallflag < 100)) {

		if (stallflag % 10 == 0) { //check for back and forthness
			for (i = 0; i < ndata; i++) {
				cycle_check_assign[i] = cluster_assign[i];
			}
		}*/

	/*	printf("\n");
		for (i = 0; i < ndata; i++) {
			printf("%d,", cluster_assign[i]);
		}
		printf("\n");
		printf("\n");
		for (i = 0; i < ndata*dim; i+=dim) {
			printf("%d) ", i/dim);
			for (j = 0; j < dim; j++) {
				printf("%f,", data[i+j]);
			}
			printf("\t%d\n", cluster_assign[i/dim]);
		}*/

		//initialize array for comparison
	/*	for (i = 0; i < ndata; i++) {
			temp_assign[i] = cluster_assign[i];
		}
		//do "next" assignment. if nothing changes, we will be stopping
		empty_idx = assignData(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign);

		while ((empty_idx < k) && (f < k)) {
			f++;
			initializeCentroids(dim, ndata, data, k, cluster_centroid, empty_idx, k, k); 
			empty_idx = assignData(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign);
		}	
		f = 0;
		
		doneflag = 1;
		for (i = 0; i < ndata; i++) {
			if (temp_assign[i] != cluster_assign[i]) {
				doneflag = 0;
				//only need to calculate new centroids if a data point doesn't match
				calculateCentroids(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign);
				break;
			}
		}

		for (i = 0; i < ndata; i++) {
			if (cluster_assign[i] != cycle_check_assign[i]) {
				break;
			}
		}

		//if it got to the end of the loop without breaking
		if (i == ndata) {
			cycle_check++;
		}

		if (cycle_check > 0) {
			doneflag = 1;
		}*/

	/*	printf("centroids\n");
		for (i = 0; i < k; i++) {
			for (j = 0; j < dim; j++) {
				printf("%f,", cluster_centroid[i][j]);
			}
			printf("\n");
		}*/
	/*	stallflag++;
		cycle_cnt++;
	}

	for (i = 0; i < k; i++) {
		loopstart = cluster_start[i]*dim;
		loopend =  loopstart + cluster_size[i]*dim;
		for (j = loopstart; j < loopend; j+=dim) {
			//get distance to each point in each cluster and remember the largest
			temp_max = calculateClusterDistance(dim, j, data, i, cluster_centroid);
			max_dist = MAX(max_dist, temp_max);
		}
		cluster_radius[i] = max_dist;
		
	}*/
/*	printf("\n");
	for (i = 0; i < ndata; i++) {
		printf("%d,", cluster_assign[i]);
	}
	printf("\n");
	printf("\n");
	for (i = 0; i < ndata*dim; i+=dim) {
		printf("%d) ", i/dim);
		for (j = 0; j < dim; j++) {
			printf("%f,", data[i+j]);
		}
		printf("\t%d\n", cluster_assign[i/dim]);
	}*/
/*	free(temp_assign);	

	return cycle_cnt;
}*/

void calculateCentroids(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign) {
	//keep in mind cluster_start values are in terms of cluster_assign
	int i, j, l, loopstart, loopend;
	double temp_centroid[DIM] = {0.0}; 
	
	for (i = 0; i < k; i++) {
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
		//printf("cluster centroid ");
		for (j = 0; j < dim; j++) {
			//printf("%f\t,", temp_centroid[j]); 
			cluster_centroid[i][j] = temp_centroid[j] / cluster_size[i];
			temp_centroid[j] = 0.0;
		}
		//printf("\n");
	}	


}

void calculateInitCentroids(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int clust1, int clust2) {
	//keep in mind cluster_start values are in terms of cluster_assign
	int i, j, l, loopstart, loopend, range;
	double temp_centroid[DIM] = {0.0}; 
	
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

double calculatePointDistance(int dim, int clust_idx, int query_idx, double *query, double **cluster_centroid) {
	int i;
	double distance = 0.0;

	for (i = 0; i < dim; i++) {
		distance += pow(cluster_centroid[clust_idx][i] - query[query_idx+i], 2);
	}

	return distance;
}

double distance(int dim, int query_idx, double *query, int data_idx, double *data) {
	int i;
	double distance = 0.0;

	for (i = 0; i < dim; i++) {
		distance += pow(data[data_idx+i] - query[query_idx+i], 2);
	}

	return distance;
}

int assignData(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign) {
	int i, j, l, cent_idx, temp_data_idx = 0, clust_size_idx_cnt = 0, empty_idx = k+5, temp_size = 0;
	double temp_distance, distance = 10000.0;
	double *temp_data, *temp_assign;

	temp_data = (double *)malloc(sizeof(double) * ndata * dim);
	temp_assign = (double *)malloc(sizeof(double) * ndata);

	//intitialize temp data
	for (i = 0; i < ndata*dim; i++) {
		temp_data[i] = data[i];
	}

	//initialize cluster size
	for (i = 0; i < k; i++) {
		cluster_size[i] = 0;
	}


	for (i = 0; i < ndata*dim; i+=dim) { //each data point
		for (j = 0; j < k; j++) { //each centroid 
			temp_distance = calculateClusterDistance(dim, i, data, j, cluster_centroid);

			distance = MIN(distance, temp_distance);

			if (distance == temp_distance)
				cent_idx = j;
		}
		distance = 10000.0;
		cluster_assign[i/dim] = cent_idx;
	}

	for (l = 0; l < k; l++) { // for each cluster 
		for (i = 0; i < ndata*dim; i+=dim) { //for each data point	
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
	for (i = 0; i < ndata; i++) {
		cluster_assign[i] = temp_assign[i];
	}

	//get cluster_starts from sizes
	//printf("cluster size, ");
	for (l = 0; l < k; l++) {
		cluster_start[l] = clust_size_idx_cnt*dim;
		clust_size_idx_cnt += cluster_size[l];
	//	printf("%d,", cluster_size[l]); 
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

int assignInitData(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int clust1, int clust2) { //ndata does not refer to whole data set, only to amount in the specified cluster
	int i, j, l, cent_idx, temp_data_idx = cluster_start[clust1], clust_size_idx_cnt = cluster_start[clust1], empty_idx = k+5, range, temp_size = 0;
	int clust1_size = 0, clust2_size = 0;
	double temp_distance, distance = 10000.0;
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
			temp_distance = calculateClusterDistance(dim, i, data, j, cluster_centroid);

			distance = MIN(distance, temp_distance);

			if (distance == temp_distance)
				cent_idx = j;
		}
		distance = 10000.0;
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
	//free(temp_data);

	return empty_idx;
}

int initializeCentroids(int dim, int ndata, double *data, int *cluster_size, int *cluster_start, double **cluster_centroid, int place_idx, int cent_idx, int cent_count, int *cluster_assign) {
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
			sums_hold = sums_hold / SIZEFACTOR; // for large data sets, don't want to overflow our array
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

	printf("\nCentroid round %d\n", place_idx);
	for (i = 0; i < cent_count; i++) {
		for (j = 0; j < dim; j++) {
			printf("%f,", cluster_centroid[i][j]);
		}
		printf("\n");
	}

	printf("\nsums idx %d\n", sums_idx);

	free(sums_arr);

	return sums_idx;
}




int main() {
	
	int i, j, kcheck = 1, rint, cycles, pointcnt, next_clust = 0;
	double r;
	double *data, **cluster_centroid, *cluster_radius, *query, *result;
	int *cluster_assign, *cluster_size,  *cluster_start;

	data = (double *)malloc(sizeof(double) * N * DIM);
	cluster_assign = (int *)malloc(sizeof(int) * N);
	cluster_size = (int *)malloc(sizeof(int) * K);
	cluster_start = (int *)malloc(sizeof(int) * K);
	cluster_centroid = (double **)malloc(sizeof(double *) * K);
	cluster_radius = (double *)malloc(sizeof(double) * K);

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


	srand(23);
	
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
		cluster_start[i] = N*DIM;
	}
	
	cluster_size[0] = N;
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

	printf("\n");
	for (i = 0; i < N*DIM; i+=DIM) {
		printf("%d) ", i/DIM);
		for (j = 0; j < DIM; j++) {
			printf("%f,", data[i+j]);
		}
		printf("\n");
	}

	printf("\nCentroid\n");
	for (i = 0; i < 1; i++) {
		for (j = 0; j < DIM; j++) {
			printf("%f,", cluster_centroid[i][j]);
		}
		printf("\n");
	}

	while (kcheck < K) {
		next_clust = initializeCentroids(DIM, N, data, cluster_size, cluster_start, cluster_centroid, kcheck, next_clust, kcheck+1, cluster_assign);
		kcheck++;

		
	}

	/*for (i = 0; i < Q * DIM; i++) {
		r = ((double)rand() / (double)(RAND_MAX)) * 100.0;
		query[i] = r; 
		printf("%f,", r);
	}*/


	//printf("dank\n");

	//cycles = kmeans(DIM, N, data, K, cluster_size, cluster_start, cluster_radius, cluster_centroid, cluster_assign);

	printf("\nNew Centroid\n");
	for (i = 0; i < K; i++) {
		for (j = 0; j < DIM; j++) {
			printf("%f,", cluster_centroid[i][j]);
		}
		printf("\n");
	}

	printf("\nCentroid sizes ");
	for (i = 0; i < K; i++) {
		printf("%d,", cluster_size[i]);
	}
	printf("\n");

	//printf("\ncycles %d\n", cycles);

	printf("\n");

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

	free(result);
	free(query);
	free(cluster_radius);
	free(cluster_centroid);
	free(cluster_start);
	free(cluster_size);
	free(data);

	return 0;
}


