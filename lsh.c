#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "lsh.h"
#include "util/dataFunctions.h"
#include "util/compFunctions.h"

void countSort(int **arr, int n, int idx, int m, int max) {
    int i, j;
    int **output, *count;

    count = (int *)malloc(sizeof(int) * max);
    output = (int **)malloc(sizeof(int *) * n);

    for (i = 0; i < n; i++) {
        output[i] = (int *)malloc(sizeof(int) * m);
    }

    for (i = 0; i < max; i++) {
        count[i] = 0;
    }

    for (i = 0; i < n; i++) {
        count[arr[i][idx]]++;
    }

    for (i = 1; i < max; i++) {
        count[i] += count[i - 1];
    }
    for (i = 0; i < max; i++) {
    }

    for (i = n-1; i >= 0; i--) {
        for (j = 0; j < m; j++) {
            output[count[arr[i][idx]] - 1][j] = arr[i][j];
        }
        count[arr[i][idx]]--;
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            arr[i][j] = output[i][j];
        }
    }

    for (i = 0; i < n; i++) {
        free(output[i]);
    }
    free(count);
}

void radsort(int **arr, int n, int m, int *max) {
    int i;

    for (i = m-2; i >= 0; i--) {
        countSort(arr, n, i, m, max[i]+1);
    }
}

int local_search(int dim, int ndata, int q0, double *data, int *cluster_size, int *cluster_start, int m, int *hash_vals, int *temp_hash, int running_cnt, double *query, double *result, double **cluster_bdry, double *cluster_radius, double **cluster_centroid) {
	int q, qhash, i, j, l, min_clust_idx = -1, min_point_idx, loopstart, loopend, count = 0;
	double min_distance = DBL_MAX, current_distance = 0.0, sum, dist2bound, dist2rad;
    double *clust_dist_arr;

    clust_dist_arr = (double *)malloc(sizeof(double) * running_cnt);

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
			/*for (i = 0; i < dim; i++) {
				result[q+i] = -1.0;
			}*/
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
			/*for (i = 0; i < dim; i++) {
				result[q+i] = data[min_point_idx+i];
			}*/
		}
        
        //now we must check other clusters to see if they are closer than our closest point
        for (i = 0; i < running_cnt; i++) {
            dist2bound = pnt2bdry(dim, cluster_bdry, query, q, i);
            dist2rad = pnt2centDistance(dim, i, q, query, cluster_centroid) - cluster_radius[i];
            clust_dist_arr[i] = MAX(dist2bound, dist2rad);
        }
        for (i = 0; i < running_cnt; i++) {
            if (i == min_clust_idx) {}//already checked this cluster
            else {
                loopstart = cluster_start[i];
                loopend = cluster_start[i] + cluster_size[i]*dim;
                if (clust_dist_arr[i] < min_distance) { //cluster is closer than closest point in current cluster
                    for (j = loopstart; j < loopend; j+=dim) { //check all points in this cluster
                        current_distance = pnt2pntDistance(dim, q, query, j, data);
                        min_distance = MIN(min_distance, current_distance);
                        if (min_distance == current_distance) {
                            min_point_idx = j;
                        }
                        count++;
                    }
                    
                }
            }
        }

        for (i = 0; i < dim; i++) {
            result[q+i] = data[min_point_idx+i];
        }

        //Resest distance and cluster index
        min_distance = DBL_MAX;
		min_clust_idx = -1;
	}
    free(clust_dist_arr);

    return count;
}

void getBoundries(int dim, int ndata, double *data, int num_clusters, int *cluster_size, int *cluster_start, double **cluster_bdry, double **cluster_centroid, double *cluster_radius) {
    int i, j, k, loopstart, loopend;
    double temp_min = DBL_MAX, temp_max = 0.0, temp_sum = 0.0;

    for (i = 0; i < num_clusters; i++) {
        loopstart = cluster_start[i];
        loopend = cluster_start[i] + cluster_size[i]*dim;
        //get cluster boundry values and centroid values
        for (j = 0; j < dim; j++) { // each dim
            for (k = loopstart + j; k < loopend; k+=dim) { // each data point
                temp_sum += data[k];
                temp_min = MIN(temp_min, data[k]);
                temp_max = MAX(temp_max, data[k]);
            }
            //assign values in respective arrays
            cluster_centroid[i][j] = temp_sum / cluster_size[i];
            cluster_bdry[i][j*2] = temp_min;
            cluster_bdry[i][j*2+1] = temp_max;
            //reset counting/comparison values
            temp_sum = 0; 
            temp_min = DBL_MAX;
            temp_max = 0.0;
        }
        temp_max = 0.0;
        //get cluster_radius values
        for (j = loopstart; j < loopend; j+=dim) {
            temp_max = MAX(temp_max, pnt2centDistance(dim, i, j, data, cluster_centroid));
        }
        cluster_radius[i] = temp_max;
        temp_max = 0.0;
    }
}

void rearrange_data(double *data, int *cluster_size, int *cluster_start, int *hash_assign, int *hash_vals, int **H, int running_cnt, int m, int ndata, int dim) {
	int i, j, k, temp_idx = 0, assign_idx = 0, hash_num;
    int *temp_assign, *clust_idxs;
	double *temp_data;

	temp_data = (double *)malloc(sizeof(double) * ndata * dim);
    temp_assign = (int *)malloc(sizeof(int) * ndata);
    clust_idxs = (int *)malloc(sizeof(int) * running_cnt);

    for (i = 0; i < running_cnt; i++) {
        clust_idxs[i] = cluster_start[i];
    }

	//for (i = 0; i < running_cnt; i++) { //for each hash value
		for (j = 0; j < ndata; j++) { //for each data point
			//if (hash_assign[j] == i) { // for each hash value if it matches
            hash_num = hash_assign[j];
				for (k = 0; k < dim; k++) { // each dimmension
					temp_data[clust_idxs[hash_num]+k] = data[j*dim+k];
				}
                temp_assign[assign_idx] = hash_num;
                assign_idx++;
                //temp_idx+=dim;
                clust_idxs[hash_num] += dim;
			//}
		}
	//}

	for (i = 0; i < ndata * dim; i++) {
		data[i] = temp_data[i];
	} 
    for (i = 0; i < ndata; i++) {
        hash_assign[i] = temp_assign[i];
    }

    free(clust_idxs);
    free(temp_assign);
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

int new_check(int **H,int *hash_vals, int m, int ndata, int *hash_assign, int *clust_count, int *max_arr) {
    int i, j, assign_idx = 0, cnt = 1, vals_idx = m, new_flag = 0;
    int *temp_hash;

    temp_hash = (int *)malloc(sizeof(int) * m);
    //sort the data set, allowing duplicates for now
    radsort(H, ndata, m+1, max_arr);

    //initial stuff for first hash
    for (i = 0; i < m; i++) {
        temp_hash[i] = H[0][i];
        hash_vals[vals_idx+i] = H[0][i];
    }
    hash_assign[H[0][m]] = assign_idx;
    clust_count[assign_idx]++;

    //get rid of duplicates and add the hashes to hash_vals
    for (i = 1; i < ndata; i++) {
        new_flag = 0;
        for (j = 0; j < m; j++) {
            if (H[i][j] != temp_hash[j]) {
                new_flag = 1;
                break;
            }
        }
        if (new_flag == 1) { //new hash
            for (j = 0; j < m; j++) {
                temp_hash[j] = H[i][j];
                hash_vals[vals_idx+j] = H[i][j];
            }
            assign_idx++;
            vals_idx+=m;
            cnt++;
            hash_assign[H[i][m]] = assign_idx;
            clust_count[assign_idx]++;
        }
        else { //duplicate hash
            hash_assign[H[i][m]] = assign_idx;
            clust_count[assign_idx]++;
        }
    }

    return cnt;
}

int LSH(int dim, int ndata, double *data, int m, double **r, double *b, double w, int num_clusters, int *cluster_size, int *cluster_start, int **H, int *hash_vals, int *hash_assign) {
	
	int i, j, k, h_i, hash_x_idx = 0, hash_y_idx = 0;
	int H_count = 0, running_cnt = 0, clust_start_idx = 0;
	int *clust_cnt, *max_arr;
	
	clust_cnt = (int *)malloc(sizeof(int) * ndata); //count of how many times that hash value appears
    max_arr = (int *)malloc(sizeof(int) * m);

	for (i = 0; i < ndata * m; i++) {
		hash_vals[i] = -1;
	}
	for (i = 0; i < ndata; i++) {
		clust_cnt[i] = 0;
	}
    for (i = 0; i < m; i++) {
        max_arr[i] = 0;
    }

	//get values for each data point for each of the m vectors
	for (i = 0; i < ndata*dim; i+=dim) { //for each data point
		for (j = 0; j < m; j++) { //for each m value
			h_i = (int) (fabs(dotprod(dim, data, r, i, j) - b[j])) / w; //because int division, we get the floor
			H[hash_x_idx][hash_y_idx] = h_i;
			hash_y_idx++;
            max_arr[j] = MAX(max_arr[j], h_i);
		}
		hash_y_idx = 0;
		hash_x_idx++;
	}

    running_cnt = new_check(H, hash_vals, m, ndata, hash_assign, clust_cnt, max_arr);

	//compare hash values to get the number of clusters
	/*for (i = 0; i < ndata; i++) {
		running_cnt = check_hash(H, hash_vals, clust_cnt, i, m, running_cnt, hash_assign);
	}*/

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
    free(max_arr);

	return running_cnt;
}

void runLSH(char *path, int ndata, int dim, int m, int w, int q, double *query, double *result) {
	
	int qi, i, j, num_clusters = 0, sum = 0, temp_hash_sum = 0, q_idx = 0, count, qcnt;
	int *cluster_size, *cluster_start, **H, *hash_vals, *query_hash, *hash_assign;
    float *ft_data;
    double rnum, cluster_time_used, search_time_used;
	double *data, **r, *b, *distance_arr, **cluster_bdry, **cluster_centroid, *cluster_radius;
    clock_t startC, endC, startS, endS;

    ft_data = (float *)malloc(sizeof(float) * ndata * dim);
	data = (double *)malloc(sizeof(double) * ndata * dim);
	b = (double *)malloc(sizeof(double) * m);
	cluster_size = (int *)malloc(sizeof(int) * ndata);
	cluster_start = (int *)malloc(sizeof(int) * ndata);
	hash_vals = (int *)malloc(sizeof(int) * m * ndata);
	query_hash = (int *)malloc(sizeof(int) * m * q);
    hash_assign = (int *)malloc(sizeof(int) * ndata);
    distance_arr = (double *)malloc(sizeof(double) * q);
	H = (int **)malloc(sizeof(int *) * ndata);
	r = (double **)malloc(sizeof(double *) * m);

	//initialize 2d array
	for (i = 0; i < ndata; i++) {
		H[i] = (int *)malloc(sizeof(int) * (m+1));
	}
	for (i = 0; i < m; i++) {
		r[i] = (double *)malloc(sizeof(double) * dim);
	}

	//initialize values of size and start
	for (i = 0; i < ndata; i++) {
		cluster_size[i] = 0;
		cluster_start[i] = -1;
	}

    for (i = 0; i < ndata; i++) {
        H[i][m] = i;
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

    printf("\nStarting LSH algorithm\n");
    startC = clock();
	num_clusters = LSH(dim, ndata, data, m, r, b, w, num_clusters, cluster_size, cluster_start, H, hash_vals, hash_assign);

    printf("\n%d\n", num_clusters);
    //after the number of clusters is known we can allocate space for boundries and radii
    cluster_bdry = (double **)malloc(sizeof(double *) * num_clusters);
    cluster_centroid = (double **)malloc(sizeof(double *) * num_clusters);
    cluster_radius = (double *)malloc(sizeof(double) * num_clusters);
    
    for (i = 0; i < num_clusters; i++) {
        cluster_bdry[i] = (double *)malloc(sizeof(double) * 2 * dim);
        cluster_centroid[i] = (double *)malloc(sizeof(double) * dim);
    }   
    
    //get the boundries and radii for exact cluster searching
    getBoundries(dim, ndata, data, num_clusters, cluster_size, cluster_start, cluster_bdry, cluster_centroid, cluster_radius);
    endC = clock();
    cluster_time_used = ((double) (endC - startC)) / CLOCKS_PER_SEC;
    printf("\n------------------\nThe Clustering time used by the algorithm was %f seconds\n------------------\n\n", cluster_time_used);
    /*for (i = 0; i < 50; i++) {
        for (j = 0; j < dim; j++) {
            printf("%f-%f, ", cluster_bdry[i][j*2], cluster_bdry[i][j*2+1]);
        }
        printf("\n");
    }*/

	//printf("\nThe number of clusters is %d\n", num_clusters);

    //write_results(dim, ndata, data, hash_assign);
	
    /*printf("Cluster Sizes:\n");
    for (i = 0; i < num_clusters; i++) {
        printf("%d\n", cluster_size[i]);
    }
    printf("\n");*/

    /*for (i = 0; i < num_clusters*m; i+=m) {
        printf("%d) ", i/m);
        for (j = 0; j < m; j++) {
            printf("%d-", hash_vals[i+j]);
        }
        printf("\t");
    }
    printf("\n");*/

    /*for (i = 0; i < 100 * m; i+=m) {
        printf("%d) <%d> ", i/m, hash_assign[i/m]);
        for (j = 0; j < m; j++) {
            printf("%d-", H[i/m][j]);
        }
        printf("\n");
    }*/

    /*for (i = 0; i < 100 * dim; i+=dim) {
        printf("%f %f %d %d\n", data[i], data[i + 1], H[i/dim][0], H[i/dim][1]);
    }*/

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
    for (qcnt = 10; qcnt < 10001; qcnt*=10) {
    startS = clock();
	count = local_search(dim, ndata, qcnt, data, cluster_size, cluster_start, m, hash_vals, query_hash, num_clusters, query, result, cluster_bdry, cluster_radius, cluster_centroid);
    endS = clock();
    search_time_used = ((double) (endS - startS)) / CLOCKS_PER_SEC;
    printf("\n------------------\nThe Searching time used by the algorithm for %d q's was %f seconds\n------------------\n\n", qcnt, search_time_used);

    }
    //printf("Number of points checked %d\n\n", count);

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
        //printf("%f\n", distance_arr[i/dim]);
    }
    //printf("\n");


	//free memory
	for (i = 0; i < ndata; i++) {
		free(H[i]);
	}
	for (i = 0; i < m; i++) {
		free(r[i]);
	}
    for (i = 0; i < num_clusters; i++) {
        free(cluster_bdry[i]);
        free(cluster_centroid[i]);
    }
    free(cluster_radius);
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
