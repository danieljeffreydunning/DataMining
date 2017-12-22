#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "lsh.h"

#define N 16000
#define W 10000
#define M 3
#define DIM 8
#define Q 1

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

//function prototypes
int LSH(int dim, int ndata, double *data, int m, double **r, double *b, double w, int num_clusters, int *cluster_size, 
	int *cluster_start, int **H, int *hash_vals);
int dotprod(int dim, double *data, double **r, int d_idx, int r_idx); //dot product
int check_hash(int **H, int *hash_vals, int *clust_cnt, int idx, int m, int running_cnt, int *hash_assign); //compares a new hash value to existi											ng ones. returns total number of clusters we have so far
void rearrange_data(double *data, int *cluster_size, int *cluster_start, int *hash_assign, int *hash_vals, int **H, int running_cnt, int m, int ndata, int dim); 
void local_search(int dim, int ndata, double *data, int *cluster_size, int *cluster_start, int m, int *hash_vals, int *temp_hash, int running_cnt, double *query, double *result);
double distance(int dim, int query_idx, double *query, int data_idx, double *data);

/*******************/

double distance(int dim, int query_idx, double *query, int data_idx, double *data) {
	int i;
	double distance = 0.0;

	for (i = 0; i < dim; i++) {
		distance += pow(data[data_idx+i] - query[query_idx+i], 2);
	}

	return distance;
}

void local_search(int dim, int ndata, double *data, int *cluster_size, int *cluster_start, int m, int *hash_vals, int *temp_hash, int running_cnt, double *query, double *result) {
	int q, i, j, l, min_clust_idx = -1, min_point_idx;
	double min_distance = 1000000.0, current_distance = 0.0, sum;

	for (q = 0; q < Q*m; q+=m) {
		for (i = 0; i < m*running_cnt; i+=m) {
			min_clust_idx = i / m;
			for (j = 0; j < m; j++) {
				if (temp_hash[q+j] != hash_vals[i+j]) {
					printf("%d <-> %d\n", temp_hash[q+j], hash_vals[i+j]);
					min_clust_idx = -1;
					break;
				}
			}
			if (min_clust_idx != -1) {
				printf("\n%d\n", min_clust_idx);
				break;
			}
		}

		//after getting closest cluster

		if (min_clust_idx == -1) { //did not match a cluster hash
			for (i = 0; i < dim; i++) {
				result[q+i] = -1.0;
			}
		}
		else { //search in matched cluster
			//calculate distance to each point in the closest cluster
			for (l = cluster_start[min_clust_idx]; l < cluster_start[min_clust_idx]+cluster_size[min_clust_idx]*dim; l+=dim) {
				current_distance = distance(dim, q, query, l, data);
				//count++;
				//printf("%f\n", current_distance);
				
				min_distance = MIN(min_distance, current_distance);
				if (min_distance == current_distance) {
					min_point_idx = l;
				}	
			}
			for (i = 0; i < dim; i++) {
				result[q+i] = data[min_point_idx+i];
			}
		}
		min_clust_idx = -1;
	}
}

void rearrange_data(double *data, int *cluster_size, int *cluster_start, int *hash_assign, int *hash_vals, int **H, int running_cnt, int m, int ndata, int dim) {
	int i, j, k, temp_idx = 0;
	double *temp_data;

	temp_data = (double *)malloc(sizeof(double) * ndata * dim);

	for (i = 0; i < running_cnt; i++) { //for each hash value
		for (j = 0; j < ndata; j++) { //for each data point
			if (hash_assign[j] == i) { // for each hash value if it matches
				for (k = 0; k < dim; k++) { // each dimmension
					temp_data[temp_idx*dim+k] = data[j*dim+k];
				}
			}
		}
	}

	for (i = 0; i < ndata * dim; i++) {
		data[i] = temp_data[i];
	} 

	free(temp_data);
}

int check_hash(int **H, int *hash_vals, int *clust_cnt, int idx, int m, int running_cnt, int *hash_assign) {
	int i, j, match_flag = 1, arr_val;
	char char_flag = 1;

	if (running_cnt == 0) { //base case, so we don't have to worry about non-indexes
		for (i = 0; i < m; i++) {
			hash_vals[running_cnt*m+i] = H[idx][i];
		}
		hash_assign[idx] = 0;
		arr_val = clust_cnt[running_cnt];
		arr_val++;
		clust_cnt[running_cnt] = arr_val;
		running_cnt++;
		return running_cnt;
	}

	//printf("%d\n", running_cnt);
	for (i = 0; i < running_cnt; i+=m) { //for each item in our hash values
		for (j = 0; j < m; j++) { //for each of m values in each hash value
			if (H[idx][j] != hash_vals[i+j]) { //if it is not a match
				match_flag = 0;
				break; 	
			}
		} 
		if (match_flag > 0) { // if we found a match
			arr_val = clust_cnt[i/m];
			arr_val++;
			clust_cnt[i/m] = arr_val;
			hash_assign[idx] = i / m;
			return running_cnt; //return, don't need to check the other hash values, running_cnt stays the same
		}
		match_flag = 1;
	}

	//if the hash did not match any others it is added to the list of hashes and the number of hashes increases by one
	for (i = 0; i < m; i++) {
		hash_vals[running_cnt*m+i] = H[idx][i];
	}
	clust_cnt[running_cnt]++;
	hash_assign[idx] = running_cnt; 
	running_cnt++;
	return running_cnt;
}

int dotprod(int dim, double *data, double **r, int d_idx, int r_idx) {
	int i, sum = 0;

	for (i = 0; i < dim; i++) {
		sum += data[d_idx+i] * r[r_idx][i];
	}

	//printf("Sum: %d\n", sum);
	return sum;
}

int LSH(int dim, int ndata, double *data, int m, double **r, double *b, double w, int num_clusters, int *cluster_size,
	int *cluster_start, int **H, int *hash_vals) {
	
	int i, j, k, h_i, hash_x_idx = 0, hash_y_idx = 0;
	int H_count = 0, running_cnt = 0, clust_start_idx = 0;
	int *clust_cnt, *tmp_size, *temp_start, *hash_assign;
	//double *temp_d, *temp_r;
	
	//hash_vals = (int *)malloc(sizeof(int) * ndata * m); //storage of each hash value
	clust_cnt = (int *)malloc(sizeof(int) * ndata); //count of how many times that hash value appears
	hash_assign = (int *)malloc(sizeof(int) * ndata); //which hash value each data point is assigned
	//temp_d = (double *)malloc(sizeof(double) * dim);
	//temp_r = (double *)malloc(sizeof(double) * dim);

	for (i = 0; i < ndata * m; i++) {
		hash_vals[i] = -1;
	}
	for (i = 0; i < ndata; i++) {
		clust_cnt[i] = 0;
	}

	//get values for each data point for each of the m vectors
	for (i = 0; i < ndata*dim; i+=dim) { //for each data point
		for (j = 0; j < m; j++) { //for each m value
			//temp_d = &data[i];
			//temp_r = &r[j][0];

			h_i = (dotprod(dim, data, r, i, j) - b[j]) / w; //because int division, we get the floor
			H[hash_x_idx][hash_y_idx] = h_i;
			hash_y_idx++;
		}
		hash_y_idx = 0;
		hash_x_idx++;
	}
	/*printf("H\n");
	for (i = 0; i < ndata; i++) {
		printf("%d) ", i);
		for (j = 0; j < m; j++) {
			printf("%d ", H[i][j]);
		}
		printf("\n");
	}
	printf("\n");*/

	//compare hash values to get the number of clusters
	for (i = 0; i < ndata; i++) {
		running_cnt = check_hash(H, hash_vals, clust_cnt, i, m, running_cnt, hash_assign);
	}

	/*printf("hash_vals\n");
	for (i = 0; i < running_cnt*m; i++) {
		printf("%d,", hash_vals[i]);
	}
	printf("\n");*/

	for (i = 0; i < running_cnt; i++) {
		cluster_size[i] = clust_cnt[i];
	}
	cluster_start[0] = 0;
	clust_start_idx += cluster_size[0];
	for (i = 1; i < running_cnt; i++) {
		cluster_start[i] = cluster_size[clust_start_idx];
		clust_start_idx += cluster_size[i];
	}

	/*printf("hash_assign\n");
	for (i = 0; i < ndata; i++) {
		printf("%d, ", hash_assign[i]);
	}
	printf("\n\n");*/

	for (i = 0; i < running_cnt; i++) {
		printf("%d\n", cluster_size[i]);
	}

	rearrange_data(data, cluster_size, cluster_start, hash_assign, hash_vals, H, running_cnt, m, ndata, dim);

	return running_cnt;
}

int main(int argc, char** argv) {
	
	int q, i, j, num_clusters = 0, w, m, ndata, dim, k, sum = 0, temp_hash_sum = 0, q_idx = 0;
	double rnum;
	int *cluster_size, *cluster_start, **H, *hash_vals, *query_hash;
	double *data, **r, *b, *query, *result;

	ndata = N;
	dim = DIM;
	m = M;
	w = W;

	data = (double *)malloc(sizeof(double) * ndata * dim);
	b = (double *)malloc(sizeof(double) * m);
	cluster_size = (int *)malloc(sizeof(int) * ndata);
	cluster_start = (int *)malloc(sizeof(int) * ndata);
	hash_vals = (int *)malloc(sizeof(int) * m * ndata);
	query_hash = (int *)malloc(sizeof(int) * m * Q);
	query = (double *)malloc(sizeof(double) * Q * dim);
	result = (double *)malloc(sizeof(double) * Q * dim);
	H = (int **)malloc(sizeof(int *) * ndata);
	r = (double **)malloc(sizeof(double *) * m);

	//initialize 2d array
	for (i = 0; i < ndata; i++) {
		H[i] = (int *)malloc(sizeof(int) * m);
	}
	for (i = 0; i < m; i++) {
		r[i] = (double *)malloc(sizeof(double) * dim);
	}

	//initialize values of size and start
	for (i = 0; i < ndata; i++) {
		cluster_size[i] = 0;
		cluster_start[i] = -1;
	}

	srand(69);

	//initialize random data points
	//printf("\nData\n");
	for (i = 0; i < ndata*dim; i+=dim) {
		//printf("%d) ", i/dim);
		for (j = 0; j < dim; j++) {
			rnum = ((double)rand() / (double)(RAND_MAX)) * 100.0;
			data[i+j] = rnum;
			//printf("%f,", rnum);
		}
		//printf("\n");
	}

	//initialize m random vectors
	printf("\nR vectors\n");
	for (i = 0; i < m; i++) {
		printf("%d) ", i);
		for (j = 0; j < dim; j++) {
			rnum = ((double)rand() / (double)(RAND_MAX)) * 100.0;
			r[i][j] = rnum;
			printf("%f,", rnum);
		}
		printf("\n");
	}
	
	//initialize random queries
	printf("\nQueries\n");
	for (i = 0; i < Q*dim; i+=dim) {
		printf("%d) ", i/dim);
		for (j = 0; j < dim; j++) {
			rnum = ((double)rand() / (double)(RAND_MAX)) * 100.0;
			query[i+j] = rnum;
			printf("%f,", rnum);
		}
		printf("\n");
	}

	for (i = 0; i < m; i++) {
		b[i] = 0;
	}


	num_clusters = LSH(dim, ndata, data, m, r, b, w, num_clusters, cluster_size, cluster_start, H, hash_vals);

	printf("The number of clusters is %d\n\n", num_clusters);
	
	for (i = 0; i < num_clusters; i++) {
		sum += cluster_size[i];
	}

	printf("Average points per cluster is %d\n\n", sum / num_clusters);

	//get hash values for the queries
	for (q = 0; q < Q*dim; q+=dim) {
		for (i = 0; i < m; i++) {
			for (j = 0; j < dim; j++) {
				temp_hash_sum += query[q+j] * r[i][j];
			}
			temp_hash_sum = (temp_hash_sum - b[i]) / w;
			query_hash[q_idx+i] = temp_hash_sum;
			temp_hash_sum = 0;
		}
		q_idx += m;
	}

	printf("\n");
	for (i = 0; i < m; i++) {
		printf("%d, ", query_hash[i]);	
	}
	printf("\n");

	local_search(dim, ndata, data, cluster_size, cluster_start, m, hash_vals, query_hash, num_clusters, query, result);

	printf("\nResults\n");
	for (i = 0; i < Q*dim; i+=dim) {
		printf("%d) ", i/dim);
		for (j = 0; j < dim; j++) {
			printf("%f,", result[i+j]);
		}
		printf("\n");
	}

	//free memory
	for (i = 0; i < ndata; i++) {
		free(H[i]);
	}
	for (i = 0; i < m; i++) {
		free(r[i]);
	}
	free(query_hash);
	free(hash_vals);
	free(query);
	free(result);
	free(cluster_size);
	free(cluster_start);
	free(b);
	free(data);

	return 0;
}
