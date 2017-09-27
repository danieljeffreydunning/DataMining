#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 5000
#define K 64 //even power of 2 
#define DIM 2
#define Q 10

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define LOG2(X) log((X)) / log(2)




int pointDistance(int dim, double *data, double *query, int comp0, int i0) {
	int i, j;
	double distance = 0.0;

	for (i = i0, j = comp0; i < i0 + dim; i++, j++) {
		distance += pow(query[j] - data[i], 2); 
	}

	distance = sqrt(distance);

	return distance;
}

/*
dim - dimensions, ndata - number of data points, data - array size dim*ndata of data (0-dim-1,...,ndata*dim-1),
k - number of cluster, cluster_size - array size k of size of each cluster, cluster_start - array size k of start idxs of each cluster,
cluster_bdry - array size 4xK (?) for cluster boundries, cluster_centroid - array size 2xK for centroid of each cluster, 
cluster_assign - array size dim*ndata to hold data points when they're being assigned to cluster after partition    
*/ 

//similar arguments for bipartitian and search

void bipartition(int dim, int i0, int im, double *data, int cluster_size[2], int cluster_start[2], double cluster_bdry[2*2*DIM], double cluster_centroid[2*2*DIM], int *cluster_assign) {
//i0 - "starting" data point, im - "ending" data point inclusive i.e for cluster size of N points, i0 = 0, im = dim * N - dim


	//double var, var_x = 0.0, var_y = 0.0, sum_x = 0, sum_y = 0, mean_x, mean_y;
	int loopstart = i0, loopend = im + 1, temp_cluster_size = im - i0 + 1, *temp_assign;	
	int i, j, k, clust1st = i0, clust2st = im, clust1sizecnt=0, clust2sizecnt=0;
	int clust1bound, clust2bound; //in array, right bound for clust1 (left is start) and left bound for clust2 (end is right) both inclusive
	//int min_y = 10000, min_x = 10000, max_y = 0, max_x = 0, minY_idx, minX_idx, maxY_idx, maxX_idx; //idx's will always be the x coordinate

	double var_arr[DIM] = {0.0}, sum_arr[DIM] = {0.0}, mean_arr[DIM], min_arr[DIM], max_arr[DIM] = {0.0};
	double var_val = 0;
	int var_idx;

	temp_assign = (int *)malloc(sizeof(int) * N);

	//initialize min_arr to 10000
	for (i = 0; i < dim; i++) {
		min_arr[i] = 10000.0;
	}

	for (i = 0; i < N; i++) {
		temp_assign[i] = cluster_assign[i];
		//printf("%d,", temp_assign[i]);
	}
	//printf("\n");

	//calculate dim sums
	for (i = 0; i < dim; i++) {
		for (j = loopstart; j < loopend; j++) { //j is index of cluster_assign
			sum_arr[i] += data[cluster_assign[j]+i];

			//also keep track of min and max values for boundries
			min_arr[i] = MIN(min_arr[i], data[cluster_assign[j]+i]);

			max_arr[i] = MAX(max_arr[i], data[cluster_assign[j]+i]);
		}
		//calculate dim means
		mean_arr[i] = sum_arr[i] / temp_cluster_size;

		//calculate dim variances
		for (j = loopstart + i; j < loopend; j+=dim) {
			var_arr[i] += pow(data[cluster_assign[j]+i] - mean_arr[i], 2);
		}
	}

	//choose variance and save index
	for (i = 0; i < dim; i++) {
		var_val = MAX(var_val, var_arr[i]);
		//save index
		if (var_val == var_arr[i])
			var_idx = i;
	}

	//partition into each side of cluster_assign based on the mean
	for (i = loopstart; i < loopend; i++) {
		//printf("i: %d\n", i);
		//printf("data[temp_assign[i]+varidx]: %f\n", data[temp_assign[i]+var_idx]);
		if (data[temp_assign[i]+var_idx] < mean_arr[var_idx]) {
			/*for (j = 0; j < dim; j++) {
				cluster_assign[clust1st+j] = data[i+j];	
			}*/
			cluster_assign[clust1st] = temp_assign[i];
			clust1bound = clust1st;
			clust1st++;
			//printf("clust1st: %d\n", clust1st);
			clust1sizecnt++;
		}
		else {
			/*for (j = 0; j < dim; j++) {
				cluster_assign[clust2st+j] = data[i+j];
			}*/
			cluster_assign[clust2st] = temp_assign[i];
			clust2bound = clust2st;
			clust2st--;
			//printf("clust2st: %d\n", clust2st);
			clust2sizecnt++;
		}
	}
	
	/*printf("\n\n");
	for (i = 0; i < temp_cluster_size; i++) {
		printf("%d,", cluster_assign[i]);
	}
	printf("\n\n");*/


	cluster_size[0] = clust1sizecnt;
	cluster_size[1] = clust2sizecnt;

	//printf("left size: %d right size: %d\n", cluster_size[0], cluster_size[1]);
	
	cluster_start[0] = i0;
	cluster_start[1] = clust2bound;

	for (j = 0; j < 2; j++) { //each cluster
		for (i = 0; i < 2*dim; i+=2) { //each dim
			if (((i/2) == var_idx) && (j == 0)) { //1st cluster variance bounds
				cluster_bdry[j*2*dim+i] = min_arr[i/2];
				cluster_bdry[j*2*dim+i+1] = mean_arr[i/2];
				
			}
			else if (((i/2) == var_idx) && (j == 1)) { //2nd cluster variance bounds
				cluster_bdry[j*2*dim+i] = mean_arr[i/2];
				cluster_bdry[j*2*dim+i+1] = max_arr[i/2];
			}
			else { //either cluster general bounds
				cluster_bdry[j*2*dim+i] = min_arr[i/2];
				cluster_bdry[j*2*dim+i+1] = max_arr[i/2];
			}
		}
	}

	/*printf("\ncluster bdry\n");
	for (i = 0; i < 2*2*dim; i++) {
		printf("%f,", cluster_bdry[i]);
	}
	printf("\n");*/

	k = 0;
	/*for (j = loopstart; j < loopend; j++) {
		data[j] = cluster_assign[k];
		k++;
	}*/

	//just for print statements
	/*printf("Mean: ");
	for (i = 0; i < dim; i++) {
		printf("%f\t", mean_arr[i]);
	}
	printf("\nVariance: ");
	for (i = 0; i < dim; i++) {
		printf("%f\t", var_arr[i]);
	}
	printf("\nStart and size: ");
	for (i = 0; i < dim; i++) {
		printf("%d,%d\t", cluster_start[i], cluster_size[i]); 
	}*/
	//printf("\nvar_idx:%d cluster1size:%d cluster2size:%d\n", var_idx, clust1sizecnt, clust2sizecnt);
	
	
	
/*NEED TO MAKE MIN AND MAX ARRAYS OF SIZE DIM (TESTING WITH 128 DIMS)*/
/*NEED TO ASSIGN BOUNDRIES AND SUCH AFTER QUESTIONS ARE ANSWERED*/


}

void biparttracker(int dim, int ndata, int depth, int cluster, int i0, int im, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_bdry, double **cluster_centroid, int *cluster_assign) {
	
	double this_cluster_boundry[2*2*DIM], this_cluster_centroid[2*2*DIM];
	int  this_cluster_size[2],  this_cluster_start[2], clust1im, clust2im, *temp_assign;
	int i, j;

	temp_assign = (int *)malloc(sizeof(int) * ndata);

	//initial biparition, all others will be in a loop until desired K is reached	
	//printf("Cluster %d\n", cluster);

	if (depth > LOG2(k)) {
		//do nothing, reached max K
		return;
	}
	else if (depth == LOG2(k)) {
		/*for (i = cluster_start[k]; i < cluster_start[k] + cluster_size[a]; i++) {
			for (j = 0; j < dim; j++) {
				//printf("\nb:%d j:%d i:%d", b, j, i);
				//printf("\ncluster_assign[i] = %d", cluster_assign[i]);
				data[b+j] = temp_data[cluster_assign[i]+j];
				//printf("\ndata[b+j] = %f", data[b+j]);
			}
			b+=dim;
		}*/
		return;
	}
	else {

		//call biparttracker
		bipartition(dim, i0, im, data, this_cluster_size, this_cluster_start, this_cluster_boundry, this_cluster_centroid, cluster_assign);

		//simple calculation to keep a lot of arithmetic out of function call
		clust1im = this_cluster_start[0] + this_cluster_size[0] - 1;	
		clust2im = this_cluster_start[1] + this_cluster_size[1] - 1;

	
		for (i = 0; i < 2; i++) {
			if (depth == LOG2(k) - 1) {
				cluster_size[cluster*2+i] = this_cluster_size[i];
				cluster_start[cluster*2+i] = this_cluster_start[i];
			}
			for (j = 0; j < dim*2; j++) {
				cluster_bdry[cluster*2+i][j] = this_cluster_boundry[2*dim*i+j];
			}
		}

		//assign the clusters their cluster IDs for this depth in cluster assign
		/*for (i = i0; i < im; i++) {
			cluster_assign[i] = cluster * 2 + cluster_assign[i];
		}*/

		biparttracker(dim, ndata, depth+1, cluster*2, i0, clust1im, data, k, cluster_size, cluster_start, cluster_bdry, cluster_centroid, cluster_assign);	
		biparttracker(dim, ndata, depth+1, cluster*2+1, this_cluster_start[1], im, data, k, cluster_size, cluster_start, cluster_bdry, cluster_centroid, cluster_assign);	
		
	}	

} 

void kdtree(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_bdry, double **cluster_centroid, int *cluster_assign) {

	double *temp_data;
	int i, j, a, b, done = 0;


	biparttracker(dim, ndata, 0, 0, 0, ndata-1, data, k, cluster_size, cluster_start, cluster_bdry, cluster_centroid, cluster_assign);
		
	//store data in temp data so we can rearrange data
	temp_data = (double *)malloc(sizeof(double) * ndata * dim);
	for (i = 0; i < ndata * dim; i++) {
		temp_data[i] = data[i];	
	}
	
	a = 0; //cluster index
	b = 0; //data array index
	//rearrange data based on cluster assign. does one cluster per while-loop iteration
	while(done < k) {		
		//i is index of cluster_assign, getting its bounds from cluster start and cluster size
		//j is each dimension of the index : b for data, and cluster_assign[i] for temp_data
		for (i = cluster_start[a]; i < cluster_start[a] + cluster_size[a]; i++) {
			for (j = 0; j < dim; j++) {
				//printf("\nb:%d j:%d i:%d", b, j, i);
				//printf("\ncluster_assign[i] = %d", cluster_assign[i]);
				data[b+j] = temp_data[cluster_assign[i]+j];
				//printf("\ndata[b+j] = %f", data[b+j]);
			}
			b+=dim;
		}
		done++;
		a++;
	}

	printf("\n");

/*	printf("cluster assign ");
	for (i = 0; i < ndata; i++) {
		printf("%d,", cluster_assign[i]);
	}
	printf("\n");

	printf("old cluster_start: ");
	for (i = 0; i < k; i++) {
		printf("%d,", cluster_start[i]); 
	}
	printf("\n"); 

	printf("cluster_size: ");
	for (i = 0; i < k; i++) {
		printf("%d,", cluster_size[i]); 
	}


	//after kdtree is complete, need to make cluster start for data array, not cluster assign
	printf("new cluster_start: ");*/
	for (i = 0, j = 0; i < k; i++) {
		cluster_start[i] = j;
		j += cluster_size[i] * dim;
	//	printf("%d,", cluster_start[i]); 
	}
	//printf("\n"); 

	free(temp_data);
}


int search_kdtree(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_bdry, double *query, double *result_pt) {

	int i, j, q, qidx, min_clust_idx, point_idx, point_check_cnt = 0;
	double q0, min0, max0, calcdist, *clust_dist_arr, compare_dist = 10000.0, temp_dist = 10000.0;

	clust_dist_arr = (double *)malloc(sizeof(double) * k);

	for (q = 0; q < Q*dim; q+=dim) { //for each query
		qidx = q;
		printf("distance array: ");
		for (i = 0; i < k; i++) { //for each cluster
			calcdist = 0;
			for (j = 0; j < 2*dim; j+=2) { //for each dimension
				q0 = query[qidx+(j/2)];
				min0 = cluster_bdry[i][j];
				//printf("min %f\n", min0);
				max0 = cluster_bdry[i][j+1];
				//printf("max %f\n", max0);
				if (q0 < min0) { //smaller than this dim min
					calcdist += pow(q0 - min0, 2);
				}
				else if (q0 > max0) { //greater than this dim max
					calcdist += pow(q0 - max0, 2);
				}
				else {} //inside min and max for this dim so "distance to cluster boundries for this dim" is 0
			}
			calcdist = sqrt(calcdist);
			clust_dist_arr[i] = calcdist;
			printf("%f,", calcdist);

			compare_dist = MIN(compare_dist, calcdist); //find smallest distance to a cluster
			if (compare_dist == calcdist) {
				//we need the index of this cluster
				min_clust_idx = i;
			}
		}
		printf("\n");
		printf("closest cluster: %d\n", min_clust_idx);
	
		//find closest point in the cluster we're in or closest to
		//compare dist now changes its meaning from distance to cluster to distance to points
		compare_dist = 10000.0;
		printf("distance: ");
		for (i = cluster_start[min_clust_idx]; i < cluster_size[min_clust_idx]*dim + cluster_start[min_clust_idx]; i+=dim) {
			temp_dist = pointDistance(dim, data, query, q*dim, i);	
			printf("%f,", temp_dist);
			compare_dist = MIN(temp_dist, compare_dist);
			if (compare_dist == temp_dist) {
				point_idx = i;
			} 			
			point_check_cnt++; //checked a new point
		}
		printf("\nclosest point in cluster: %d\n", point_idx / dim);
		//printf("i: %d\n", i);

		//check if any clusters are closer than the given point
		for (i = 0; i < k; i++) {
			if (i == min_clust_idx) {}//already checked this cluster
			else {
				if (clust_dist_arr[i] < compare_dist) { //cluster is closer than closest point in current cluster
					for (j = cluster_start[min_clust_idx]; j < cluster_size[min_clust_idx]*dim + cluster_start[min_clust_idx]; j+=dim) { //check all points in this cluster
						temp_dist = pointDistance(dim, data, query, q*dim, j);	
						//printf("%f,", temp_dist);
						compare_dist = MIN(temp_dist, compare_dist);
						if (compare_dist == temp_dist) {
							point_idx = j;
						} 			
						point_check_cnt++;
					}
					
				}
			}
		}

		printf("[new] closest point in cluster: %d\n\n", point_idx / dim);

		for (i = 0; i < dim; i++) {
			result_pt[qidx+i] = data[point_idx+i];
		}
		
	}

	

	return point_check_cnt / Q;
}

int main() {
	
	int i, j, avg_checks;
	double r;
	double *data, **cluster_boundry, **cluster_centroid, *query, *result;
	int *cluster_assign, *cluster_size,  *cluster_start;

	data = (double *)malloc(sizeof(double) * N * DIM);
	cluster_assign = (int *)malloc(sizeof(int) * N);
	cluster_size = (int *)malloc(sizeof(int) * K);
	cluster_start = (int *)malloc(sizeof(int) * K);
	cluster_boundry = (double **)malloc(sizeof(double *) * K);
	cluster_centroid = (double **)malloc(sizeof(double *) * K);

	query = (double *)malloc(sizeof(double) * Q * DIM);
	result = (double *)malloc(sizeof(double) * DIM);

	//initialize 2d arrays
	for (i = 0; i < K; i++) {
		cluster_boundry[i] = (double *)malloc(sizeof(double) * DIM * 2);
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

	for (i = 0; i < Q * DIM; i++) {
		r = ((double)rand() / (double)(RAND_MAX)) * 100.0;
		query[i] = r; 
		//printf("%f,", r);
	}

	//bipartition(DIM, 0, N*DIM - DIM, data, cluster_size, cluster_start, cluster_boundry, cluster_centroid);
	kdtree(DIM, N, data, K, cluster_size, cluster_start, cluster_boundry, cluster_centroid, cluster_assign);


	/*for (i = 0; i < K; i++) {
		for (j = 0; j < DIM*2; j+=2) {
			printf("Cluster %d Dim %d Min %f Max %f\n", i, j/2, cluster_boundry[i][j], cluster_boundry[i][j+1]);
		}
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

	printf("\nQuery: ");
	for (i = 0; i < Q*DIM; i++) {
		printf("%f,", query[i]);	
	} 
	printf("\n\n");


	avg_checks = search_kdtree(DIM, N, data, K, cluster_size, cluster_start, cluster_boundry, query, result);

	printf("Average number of checks was %d\n", avg_checks);
	printf("\n");

	return 0;
}
