#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
#include <float.h>
#include "kmeans.h"
#include "bisecting_kmeans.h"
#include "util/dataFunctions.h"
#include "util/compFunctions.h"

int search_kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, double *query, double *result_pt, int world_size, int world_rank, int q) {
	int i, j, l, min_clust_idx, min_point_idx, global_min_clust_idx, global_min_point_idx, count = 0, loopstart, loopend, newlowflag = 0, oglow = 0;
	double cent_dist = DBL_MAX, point_dist = DBL_MAX, global_cent_dist, global_point_dist, temp_dist, og_dist, rad_comp_dist;
	double *cent_comp_arr, *global_cent_comp_arr;

	cent_comp_arr = (double *)malloc(sizeof(double) * k);
    global_cent_comp_arr = (double *)malloc(sizeof(double) * k);

	for (j = 0; j < q * dim; j+=dim) { // for each query
		//initialize comp arr
		for (i = 0; i < k; i++) {
			cent_comp_arr[i] = 0.0;
		}
		for (i = 0; i < k; i++) { //for each cluster
			//temp_dist = calculatePointDistance(dim, i, j, query, cluster_centroid); //calculate distance to each centroid and get local minimum
            temp_dist = pnt2centDistance(dim, i, j, query, cluster_centroid); 

			cent_comp_arr[i] = temp_dist; //store because we'll need it later for comparing distance to radii
			
			cent_dist = MIN(cent_dist, temp_dist);
			if (cent_dist == temp_dist) {
				min_clust_idx = i;
			}
		}

		//calculate distance to each point in the closest cluster and get the local minimum
        loopstart = cluster_start[min_clust_idx];
        loopend = cluster_start[min_clust_idx] + cluster_size[min_clust_idx]*dim;
		for (l = loopstart; l < loopend; l+=dim) {
            temp_dist = pnt2pntDistance(dim, j, query, l, data);
			count++;
			point_dist = MIN(point_dist, temp_dist);
			if (point_dist == temp_dist) {
				min_point_idx = l;
			}	
		}

        //all procs need to have the same minimum distance to a point in the cluster
        MPI_Allreduce(&point_dist, &global_point_dist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        if (point_dist == global_point_dist) {
            oglow = 1;
        }
        point_dist = global_point_dist;
        og_dist = point_dist;

		//check if any cluster radius is closer than the min distance. If so, check points in that cluster
		for (i = 0; i < k; i++) {
			if (i != min_clust_idx) { //maker sure its not the closest cluster we already did
				rad_comp_dist = cent_comp_arr[i] - cluster_radius[i]; //distance to the edge of that cluster
				if (rad_comp_dist <= point_dist) {
					//calculate distance to each point in the next closest cluster(s)
					for (l = cluster_start[i]; l < cluster_start[i]+cluster_size[i]*dim; l+=dim) {
						//temp_dist = distance(dim, j, query, l, data);
                        temp_dist = pnt2pntDistance(dim, j, query, l, data);
						count++;

						point_dist = MIN(point_dist, temp_dist);
						if (point_dist == temp_dist) {
							min_point_idx = l;
						}	
					}
	
				}
			}
		}

        //if we got a new lowest, we are the new owners of the lowest
        if (point_dist < og_dist) {
            newlowflag = 1;
        }

        //all procs need to have the same minimum distance to a point in the cluster
        MPI_Allreduce(&point_dist, &global_point_dist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        //if there was a change from the original smallest
        if (og_dist != global_point_dist) {
            //if we had the shortest distance, we will add our result to the result array and send it to proc 0
            if ((newlowflag == 1) && (point_dist == global_point_dist)) { 
		        for (i = 0; i < dim; i++) {
		            result_pt[j+i] = data[min_point_idx+i];
	    	    }
                //send it to 0 if we are not 0
                if (world_rank != 0) {
                    MPI_Send(&result_pt[j], dim, MPI_DOUBLE, 0, 23, MPI_COMM_WORLD);
                }
            }
            else { //we are not owners of the shortest distance
                if (world_rank == 0) {
                    MPI_Recv(&result_pt[j], dim, MPI_DOUBLE, MPI_ANY_SOURCE, 23, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        //if there was no change
        else {
            if (oglow == 1) {
                for (i = 0; i < dim; i++) {
                    result_pt[j+i] = data[min_point_idx+i];
                }
                //send it to 0 if we are not 0
                if (world_rank != 0) {
                    MPI_Send(&result_pt[j], dim, MPI_DOUBLE, 0, 23, MPI_COMM_WORLD);
                }
            }
            else { //we are not owners of the shortest distance
                if (world_rank == 0) {
                    MPI_Recv(&result_pt[j], dim, MPI_DOUBLE, MPI_ANY_SOURCE, 23, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        oglow = 0;
        newlowflag = 0;
        point_dist = DBL_MAX;
        cent_dist = DBL_MAX;
	}

    free(global_cent_comp_arr);
    free(cent_comp_arr);

	return count;
}

int kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank) {
	int i, j, doneflag = 0, cycle_cnt = 1, loopstart, loopend, stallflag = 0, empty_idx = k+5, cycle_check = 0, f = 0, diff = 0, proc_diff = 0;
	int *temp_assign, *temp_clust_size, *cycle_check_assign, max_dist = 0.0, temp_max;
	double *temp_radius;

	temp_assign = (int *)malloc(sizeof(int) * ndata);
    temp_clust_size = (int *)malloc(sizeof(int) * k);
	cycle_check_assign = (int *)malloc(sizeof(int) * ndata);
	temp_radius = (double *)malloc(sizeof(double) * k);

	//first call of assignData
	assignData(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign, world_size, world_rank);

	MPI_Allreduce(&cluster_size[0], &temp_clust_size[0], k, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
    
	//see if there are any empty clusters
	for (i = 0; i < k; i++) {
		if (temp_clust_size[i] == 0) {
			empty_idx = i;
			break;
		}	
	}

	//make sure there are no empty clusters
	while ((empty_idx < k) && (f < k)) {
		f++;
		initializeCentroids(dim, ndata, data, k, cluster_centroid, empty_idx, k, world_rank, world_size); 
		assignData(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign, world_size, world_rank);
		empty_idx = k+5;

        //always reduce so that global sizes are stored temporarily and local sizes which are actually indexed are saved in cluster_size
	    MPI_Allreduce(&cluster_size[0], &temp_clust_size[0], k, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 

		//see if there are any empty clusters
		for (i = 0; i < k; i++) {
			if (temp_clust_size[i] == 0) {
				empty_idx = i;
				break;
			}	
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}	
	f = 0;

    //if we have no empty clusters we can calculate the new centroids
	calculateCentroids(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign, world_size, world_rank);

    //if (world_rank == 0) 
        //printf("Adjusting centroids\n...\n");

    //do the process until none of the data points change cluster (has a hardcoded limit of 50, can be modified or removed)
	while (!doneflag) {

        /*if (world_rank == 0) {
            printf("%d\n", cycle_cnt);
        }*/ 

		MPI_Barrier(MPI_COMM_WORLD);

		/*if (stallflag % 10 == 0) { //check for back and forthness
			for (i = 0; i < ndata; i++) {
				cycle_check_assign[i] = cluster_assign[i];
			}
		}*/

		//initialize array for comparison
		for (i = 0; i < ndata; i++) {
			temp_assign[i] = cluster_assign[i];
		}
		//do "next" assignment. if nothing changes, we will be stopping
		assignData(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign, world_size, world_rank);
		empty_idx = k+5;

	    MPI_Allreduce(&cluster_size[0], &temp_clust_size[0], k, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 

		//see if there are any empty clusters
		for (i = 0; i < k; i++) {
			if (temp_clust_size[i] == 0) {
				empty_idx = i;
				break;
			}	
		}

		while ((empty_idx < k) && (f < k)) {
			f++;
			initializeCentroids(dim, ndata, data, k, cluster_centroid, empty_idx, k, world_rank, world_size); 
			assignData(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign, world_size, world_rank);
			empty_idx = k+5;
            
	        MPI_Allreduce(&cluster_size[0], &temp_clust_size[0], k, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 

			//see if there are any empty clusters
			for (i = 0; i < k; i++) {
				if (temp_clust_size[i] == 0) {
					empty_idx = i;
					break;
				}	
			}

			MPI_Barrier(MPI_COMM_WORLD);
		}	
		f = 0;
		
		doneflag = 1;
		for (i = 0; i < ndata; i++) {
			if (temp_assign[i] != cluster_assign[i]) {
				diff = 1;
				//only need to calculate new centroids if a data point doesn't match
				break;
			}
		}

		MPI_Allreduce(&diff, &proc_diff, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		if (proc_diff > 0) {
			doneflag = 0;
			calculateCentroids(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign, world_size, world_rank);
		} 

        /*the following comment is to check for cycles in assignment i.e. data points switching back and forth cluster assignments*/

		/*for (i = 0; i < ndata; i++) {
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

		diff = 0;
		proc_diff = 0;
		//stallflag++;
		cycle_cnt++;
	}

    /*if (world_rank == 0) { //printf global cluster sizes
        printf("\nCluster sizes: ");
        for (i = 0; i < k; i++) {
            printf("%d ", temp_clust_size[i]);
        }
        printf("\n");
    }*/

	for (i = 0; i < k; i++) {
		loopstart = cluster_start[i];
		loopend = cluster_start[i] + cluster_size[i]*dim;
		for (j = loopstart; j < loopend; j+=dim) {
			//get distance to each point in each cluster and remember the largest
			//temp_max = calculateClusterDistance(dim, j, data, i, cluster_centroid);
            temp_max =  pnt2centDistance(dim, i, j, data, cluster_centroid);
			max_dist = MAX(max_dist, temp_max);
		}
		cluster_radius[i] = max_dist;

		max_dist = 0.0;
	}
	//reduce MAX operation for cluster radiuses
	MPI_Allreduce(&cluster_radius[0], &temp_radius[0], k, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	for (i = 0; i < k; i++) {
		cluster_radius[i] = temp_radius[i];
	}

	free(temp_assign);	
    free(temp_clust_size);
	free(cycle_check_assign);
	free(temp_radius);

	return cycle_cnt;
}

void calculateCentroids(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank) {
	//keep in mind cluster_start values are in terms of cluster_assign
	int i, j, l, loopstart, loopend;
	double *temp_centroid, *temp_proc_cent;
	int *temp_clust_size; //actually the global cluster size

    temp_centroid = (double *)malloc(sizeof(double) * dim);
    temp_proc_cent = (double *)malloc(sizeof(double) * dim);
	temp_clust_size = (int *)malloc(sizeof(int) * k); 

    for (i = 0; i < dim; i++) {
        temp_centroid[i] = 0.0;
        temp_proc_cent[i] = 0.0;
    }

    //want our temp cluster size to be the global sizes
	MPI_Allreduce(&cluster_size[0], &temp_clust_size[0], k, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
	
	for (i = 0; i < k; i++) {
		loopstart = cluster_start[i];
		loopend = cluster_start[i] + cluster_size[i]*dim;
		//for each centroid, getting the sum of all points, each of their dimensions, in temp_centroid
		for (j = loopstart; j < loopend; j+=dim) {
			for (l = 0; l < dim; l++) {	
				temp_centroid[l] += data[j+l];
			}
		} 

		MPI_Allreduce(&temp_centroid[0], &temp_proc_cent[0], dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);

		//now get the averages for each dimension in temp_centroid and put it in cluster_centroid
		for (j = 0; j < dim; j++) {
			//printf("%f\t,", temp_centroid[j]); 
			cluster_centroid[i][j] = temp_proc_cent[j] / temp_clust_size[i];
			temp_centroid[j] = 0.0;
		}
		//printf("\n");
	}	

	free(temp_clust_size);
    free(temp_proc_cent);
    free(temp_centroid);
}

void assignData(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank) {
	int i, j, l, cent_idx, temp_data_idx = 0, clust_size_idx_cnt = 0, empty_idx = k+5;
    long int l_dist, l_td;
	double temp_distance, distance = DBL_MAX;
	double *temp_data;
	int *temp_clust_size, *temp_assign;

	temp_data = (double *)malloc(sizeof(double) * ndata * dim);
	temp_assign = (int *)malloc(sizeof(int) * ndata);
	temp_clust_size = (int *)malloc(sizeof(int) * k);

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
            temp_distance = pnt2centDistance(dim, j, i, data, cluster_centroid);
            l_td = (long int) temp_distance;
            
			distance = MIN(distance, temp_distance);
            l_dist = (long int) distance;
            //use long int values to compare to avoid complications with floats
			if (l_dist == l_td)
				cent_idx = j;
		}
		cluster_assign[i/dim] = cent_idx;

        distance = DBL_MAX; //reset distance comparison
	}

	for (l = 0; l < k; l++) { // for each cluster 
		for (i = 0; i < ndata*dim; i+=dim) { //for each data point	
			//rearrange data set
			if (cluster_assign[i/dim] == l) {
				cluster_size[l]++; //keep track of each cluster size
				for (j = 0; j < dim; j++) {
					data[temp_data_idx+j] = temp_data[i+j];
				}	
				temp_assign[temp_data_idx/dim] = l;
				temp_data_idx+=dim;
			}
		}
	}

	//rearrange cluster_assign
	for (i = 0; i < ndata; i++) {
		cluster_assign[i] = temp_assign[i];
	}

	//get cluster_starts from sizes
	for (l = 0; l < k; l++) {
		cluster_start[l] = clust_size_idx_cnt;
		clust_size_idx_cnt += cluster_size[l]*dim;
	}

	free(temp_clust_size);
	free(temp_assign);
	free(temp_data);
}

void initializeCentroids(int dim, int ndata, double *data, int k, double **cluster_centroid, int place_idx, int m, int world_rank, int world_size) {
	int i, j, l, centroid_point;
    long int l_tot, l_glo, l_pnt;
	double distance = 0.0, point_min = DBL_MAX, total_max = 0.0, global_max;
	//start_idx will be 0 in initialization case, but will be the cluster to be redone if calling this function for an empty cluster

	double *temp_point;

	temp_point = (double *)malloc(sizeof(double) * dim);

	for (j = 0; j < ndata*dim; j+=dim) { //for each data point
		for (i = 0; i < m; i++) { //for given m centroids
			if (i != place_idx) {
                distance = pnt2centDistance(dim, i, j, data, cluster_centroid);
				//get min distance to a centroid for that point
				point_min = MIN(point_min, distance);
                l_pnt = (long int) point_min;
			 }
		}
		//get max min distance for each data point
		total_max = MAX(total_max, point_min);
        l_tot = (long int) total_max;
	
		//if we own the current max min, this will currently be the centroid point
		if (l_tot == l_pnt) {
			centroid_point = j;	
		}

		//reset point_min
		point_min = DBL_MAX;

	}

	//save our temp point. temp point will be used after broadcast to store value of global centroid
	for (i = 0; i < dim; i++) {
		temp_point[i] = data[centroid_point+i];
	}

	MPI_Allreduce(&total_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    l_glo = (long int) global_max;
	//if we we're the owner of global max, we had the greatest and must broadcast our data point
	if (l_glo == l_tot) {
		for (i = 0; i < world_size; i++) {
			if (i != world_rank)
				MPI_Send(&temp_point[0], dim, MPI_DOUBLE, i, world_rank, MPI_COMM_WORLD);
		}
	} 
	else {
		MPI_Recv(&temp_point[0], dim, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
		
	//now put the new centroid in its place
	for (i = 0; i < dim; i++) {
		cluster_centroid[place_idx][i] = temp_point[i];
	}

    free(temp_point);
}

void runKMeans(char *path, int ndata, int dim, int k, int q, double *query, double *result) {
	int kcheck = 1, rint, cycles, pointcnt, strsize, i, j, world_size, world_rank, proc_chunk_size, data_remainder, qcnt;
	int *cluster_assign, *cluster_size,  *cluster_start;
    float *ft_data;
    double cluster_time_used, search_time_used;
	double *proc_data, **cluster_centroid, *cluster_radius, *distance_arr;
    clock_t startC, endC, startS, endS;
    
	cluster_size = (int *)malloc(sizeof(int) * k);
	cluster_start = (int *)malloc(sizeof(int) * k);
	cluster_centroid = (double **)malloc(sizeof(double *) * k);
	cluster_radius = (double *)malloc(sizeof(double) * k);
    distance_arr = (double *)malloc(sizeof(double) * q);

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	proc_chunk_size = ndata / world_size; //the size each proc will deal with
    data_remainder = ndata % world_size;

    //"leftover" data needs to be given to procs
    if (world_rank < data_remainder)
        proc_chunk_size++;

    ft_data = (float *)malloc(sizeof(float) * proc_chunk_size * dim);
	cluster_assign = (int *)malloc(sizeof(int) * proc_chunk_size);
	proc_data = (double *)malloc(sizeof(double) * proc_chunk_size * dim);

	//initialize 2d arrays
	for (i = 0; i < k; i++) {
		cluster_centroid[i] = (double *)malloc(sizeof(double) * dim);
	}

	for (i = 0; i < k; i++) {
		cluster_size[i] = 0;
	}

    //start reading binary file
    readFloatBin(path, ft_data, proc_chunk_size * dim, world_rank);
    //convert float values into doubles because computation numbers can get quite large
    for (i = 0; i < proc_chunk_size * dim; i++) {
        proc_data[i] = (double) ft_data[i];
    }

	if (world_rank == 0) {
		//initialize centroids
		rint = rand() % proc_chunk_size;
		for (i = 0; i < dim; i++) {
			cluster_centroid[0][i] = proc_data[rint*dim+i];
		}
	
		MPI_Bcast(&cluster_centroid[0][0], dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //printf("Now initializing centroids\n...\n");
	}

	else { // not root
		MPI_Bcast(&cluster_centroid[0][0], dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

    //before while loops good to make sure everyone is on the same page
	MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank == 0) {
        printf("Starting Kmeans Algorithm\n\n");
    }

    startC = clock();
	/*while (kcheck < k) {
        //printf("%d kcheck %d\n", world_rank, kcheck);
		initializeCentroids(dim, proc_chunk_size, proc_data, k, cluster_centroid, kcheck, kcheck, world_rank, world_size);
		kcheck++;
		MPI_Barrier(MPI_COMM_WORLD);
	}*/

    runBKmeans(dim, k, proc_chunk_size, proc_data, cluster_centroid);
    //printf("here\n");
    /*for (i = 0; i < 100*dim; i+=dim) {
        for (j = 0; j < dim; j++) {
            printf("%f, ", proc_data[i+j]);
        }
        printf("\n");
    }*/
	cycles = kmeans(dim, proc_chunk_size, proc_data, k, cluster_size, cluster_start, cluster_radius, cluster_centroid, cluster_assign, world_size, world_rank);

    //printf("here\n");

    endC = clock();
    cluster_time_used = ((double) (endC - startC)) / CLOCKS_PER_SEC;
    printf("\n------------------\nThe Clustering time used by the algorithm was %f seconds\n------------------\n\n", cluster_time_used);
//if (world_rank == 0) {

	//printf("\ncycles %d\n", cycles);

    //will only work properly currently with 1 proc execution
    //write_results(dim, ndata, proc_data, cluster_assign);
//}
    for (qcnt = 10; qcnt < 10001; qcnt*=10) {
        startS = clock();
	    pointcnt = search_kmeans(dim, ndata, proc_data, k, cluster_size, cluster_start, cluster_radius, cluster_centroid, query, result, world_size, world_rank, qcnt);
        endS = clock();
        search_time_used = ((double) (endS - startS)) / CLOCKS_PER_SEC;
        printf("\n------------------\nThe Searching time used by the algorithm for %d points was %f seconds\n------------------\n\n", qcnt, search_time_used);
    }

	//printf("number of points checked %d\n\n", pointcnt);

    if (world_rank == 0) {
        for (i = 0; i < q*dim; i+=dim) {
            distance_arr[i/dim] = pnt2pntDistance(dim, i, query, i, result);
            //printf("%f\n", distance_arr[i/dim]);
        }
        //printf("\n");
    }

    for (i = 0; i < k; i++) {
        free(cluster_centroid[i]);
    }

    free(ft_data);
    free(distance_arr);
	free(result);
	free(cluster_radius);
	free(cluster_start);
	free(cluster_size);
    free(proc_data);

	MPI_Finalize();
}

