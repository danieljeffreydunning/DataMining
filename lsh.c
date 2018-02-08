#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "lsh.h"
#include "util/dataFunctions.h"
#include "util/compFunctions.h"

int local_search(int dim, int ndata, int q0, double *data, int *cluster_size, int *cluster_start, int m, int *hash_vals, int *temp_hash, int running_cnt, double *query, double *result) {
	int q, qhash, i, j, l, min_clust_idx = -1, min_point_idx, loopstart, loopend, count = 0;
	double min_distance = DBL_MAX, current_distance = 0.0, sum;

	for (q = 0, qhash = 0; q < q0*dim; q+=dim, qhash+=m) {
		for (i = 0; i < m*running_cnt; i+=m) {
			min_clust_idx = i / m;
			for (j = 0; j < m; j++) {
                //printf("%d vs %d\t", temp_hash[qhash+j], hash_vals[i+j]);
				if (temp_hash[qhash+j] != hash_vals[i+j]) {
					min_clust_idx = -1;
					break;
				}
			}
			if (min_clust_idx != -1) {
				break;
			}
		}
        /*for (i = 0; i < m; i++) {
            printf("%d-", hash_vals[min_clust_idx*m+i]);
        }
        printf("\n");*/

		//after getting closest cluster
        loopstart =  cluster_start[min_clust_idx];
        loopend = cluster_start[min_clust_idx]+cluster_size[min_clust_idx]*dim;

		if (min_clust_idx == -1) { //did not match a cluster hash
			for (i = 0; i < dim; i++) {
				result[q+i] = -1.0;
			}
            return 0;
		}
		else { //search in matched cluster
			//calculate distance to each point in the closest cluster
			for (l = loopstart; l < loopend; l+=dim) {
				current_distance = pnt2pntDistance(dim, q, query, l, data);
				min_distance = MIN(min_distance, current_distance);
				if (min_distance == current_distance) {
					min_point_idx = l;
				}
                count++;    
			}
			for (i = 0; i < dim; i++) {
				result[q+i] = data[min_point_idx+i];
			}
		}
        min_distance = DBL_MAX;
		min_clust_idx = -1;
	}
    return count;
}

void rearrange_data(double *data, int *cluster_size, int *cluster_start, int *hash_assign, int *hash_vals, int **H, int running_cnt, int m, int ndata, int dim) {
	int i, j, k, temp_idx = 0;
	double *temp_data;

	temp_data = (double *)malloc(sizeof(double) * ndata * dim);

	for (i = 0; i < running_cnt; i++) { //for each hash value
		for (j = 0; j < ndata; j++) { //for each data point
			if (hash_assign[j] == i) { // for each hash value if it matches
				for (k = 0; k < dim; k++) { // each dimmension
					temp_data[temp_idx+k] = data[j*dim+k];
				}
                temp_idx+=dim;
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
        clust_cnt[running_cnt]++;
		running_cnt++;
		return running_cnt;
	}

	for (i = 0; i < running_cnt*m; i+=m) { //for each item in our hash values
		for (j = 0; j < m; j++) { //for each of m values in each hash value
			if (H[idx][j] != hash_vals[i+j]) { //if it is not a match
				match_flag = 0;
				break; 	
			}
		} 
		if (match_flag > 0) { // if we found a match
            clust_cnt[i/m]++;
			hash_assign[idx] = i / m;
			return running_cnt; //return, don't need to check the other hash values, running_cnt stays the same
		}
		match_flag = 1; //reset flag
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

int LSH(int dim, int ndata, double *data, int m, double **r, double *b, double w, int num_clusters, int *cluster_size, int *cluster_start, int **H, int *hash_vals, int *hash_assign) {
	
	int i, j, k, h_i, hash_x_idx = 0, hash_y_idx = 0;
	int H_count = 0, running_cnt = 0, clust_start_idx = 0;
	int *clust_cnt, *tmp_size, *temp_start;
	
	clust_cnt = (int *)malloc(sizeof(int) * ndata); //count of how many times that hash value appears

	for (i = 0; i < ndata * m; i++) {
		hash_vals[i] = -1;
	}
	for (i = 0; i < ndata; i++) {
		clust_cnt[i] = 0;
	}

	//get values for each data point for each of the m vectors
	for (i = 0; i < ndata*dim; i+=dim) { //for each data point
		for (j = 0; j < m; j++) { //for each m value

			h_i = (fabs(dotprod(dim, data, r, i, j) - b[j])) / w; //because int division, we get the floor
			H[hash_x_idx][hash_y_idx] = h_i;
			hash_y_idx++;
		}
		hash_y_idx = 0;
		hash_x_idx++;
	}

	//compare hash values to get the number of clusters
	for (i = 0; i < ndata; i++) {
		running_cnt = check_hash(H, hash_vals, clust_cnt, i, m, running_cnt, hash_assign);
	}

	for (i = 0; i < running_cnt; i++) {
		cluster_size[i] = clust_cnt[i];
	}
	cluster_start[0] = 0;
	clust_start_idx += cluster_size[0]*dim;
	for (i = 1; i < running_cnt; i++) {
		cluster_start[i] = clust_start_idx;
		clust_start_idx += cluster_size[i]*dim;
	}

	/*for (i = 0; i < running_cnt; i++) {
		printf("%d\n", cluster_size[i]);
	}*/

	rearrange_data(data, cluster_size, cluster_start, hash_assign, hash_vals, H, running_cnt, m, ndata, dim);

    free(clust_cnt);

	return running_cnt;
}

void runLSH(char *path, int ndata, int dim, int m, int w, int q, double *query) {
	
	int qi, i, j, num_clusters = 0, sum = 0, temp_hash_sum = 0, q_idx = 0, count;
	int *cluster_size, *cluster_start, **H, *hash_vals, *query_hash, *hash_assign;
    float *ft_data;
    double rnum;
	double *data, **r, *b, *result, *distance_arr;

    ft_data = (float *)malloc(sizeof(float) * ndata * dim);
	data = (double *)malloc(sizeof(double) * ndata * dim);
	b = (double *)malloc(sizeof(double) * m);
	cluster_size = (int *)malloc(sizeof(int) * ndata);
	cluster_start = (int *)malloc(sizeof(int) * ndata);
	hash_vals = (int *)malloc(sizeof(int) * m * ndata);
	query_hash = (int *)malloc(sizeof(int) * m * q);
    hash_assign = (int *)malloc(sizeof(int) * ndata);
	result = (double *)malloc(sizeof(double) * q * dim);
    distance_arr = (double *)malloc(sizeof(double) * q);
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

    //start reading binary file
    readFloatBin(path, ft_data, ndata * dim, 0);
    //convert float values into doubles because computation numbers can get quite large
    for (i = 0; i < ndata * dim; i++) {
        data[i] = (double) ft_data[i];
    }

    srand(405);
	//initialize m random vectors
	for (i = 0; i < m; i++) {
		for (j = 0; j < dim; j++) {
			rnum = ((double)rand() / (double)(RAND_MAX)) * 100.0;
			r[i][j] = rnum;
		}
	}

	for (i = 0; i < m; i++) {
		b[i] = ((double)rand() / (double)(RAND_MAX)) * 5.0;;
	}

	num_clusters = LSH(dim, ndata, data, m, r, b, w, num_clusters, cluster_size, cluster_start, H, hash_vals, hash_assign);

	printf("\nThe number of clusters is %d\n", num_clusters);

    write_results(dim, ndata, data, hash_assign);
	
    printf("Cluster Sizes:\n");
    for (i = 0; i < num_clusters; i++) {
        printf("%d\n", cluster_size[i]);
    }
    printf("\n");

    /*for (i = 0; i < num_clusters*m; i+=m) {
        printf("%d) ", i/m);
        for (j = 0; j < m; j++) {
            printf("%d-", hash_vals[i+j]);
        }
        printf("\n");
    }
    printf("\n");*/

	//get hash values for the queries
	for (qi = 0; qi < q*dim; qi+=dim) {
		for (i = 0; i < m; i++) {
			for (j = 0; j < dim; j++) {
				temp_hash_sum += query[qi+j] * r[i][j];
			}
			temp_hash_sum = (temp_hash_sum - b[i]) / w;
			query_hash[q_idx+i] = temp_hash_sum;
			temp_hash_sum = 0;
		}
		q_idx += m;
	}

    /*for (i = 0; i < q*m; i+=m) {
        for (j = 0; j < m; j++) {
            printf("%d-", query_hash[i+j]);
        }
        printf("\n");
    }*/

	count = local_search(dim, ndata, q, data, cluster_size, cluster_start, m, hash_vals, query_hash, num_clusters, query, result);
    printf("Number of points checked %d\n\n", count);

    for (i = 0; i < q*dim; i+=dim) {
        /*for (j = 0; j < dim; j++) {
            printf("%f, ", result[i+j]);
        }
        printf("\n");
        for (j = 0; j < dim; j++) {
            printf("%f, ", query[i+j]);
        }
        printf("\n");*/
        distance_arr[i/dim] = pnt2pntDistance(dim, i, query, i, result);
        printf("%f\n", distance_arr[i/dim]);
    }
    printf("\n");

	//free memory
	for (i = 0; i < ndata; i++) {
		free(H[i]);
	}
	for (i = 0; i < m; i++) {
		free(r[i]);
	}
	free(query_hash);
	free(hash_vals);
    free(hash_assign);
	free(result);
    free(distance_arr);
	free(cluster_size);
	free(cluster_start);
	free(b);
    free(ft_data);
	free(data);
}
