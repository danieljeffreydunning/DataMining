#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
#include <float.h>
#include "kmeans.h"
#include "util/dataFunctions.h"
#include "util/compFunctions.h"

//function prototypes
int kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank);
void calculateCentroids(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank);
void assignData(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank);
int search_kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, double *query, double *result_pt, int world_size, int world_rank, int q);
//double calculateClusterDistance(int dim, int data_idx, double *data, int clust_idx, double **cluster_centroid);
//double calculatePointDistance(int dim, int clust_idx, int query_idx, double *query, double **cluster_centroid);
//double distance(int dim, int query_idx, double *query, int data_idx, double *data);
void initializeCentroids(int dim, int ndata, double *data, int k, double **cluster_centroid, int start_idx, int m, int world_rank, int world_size);
//end prototypes

int search_kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, double *query, double *result_pt, int world_size, int world_rank, int q) {
	int i, j, l, min_clust_idx, min_point_idx, global_min_clust_idx, global_min_point_idx, count = 0, loopstart, loopend;
	double cent_dist = DBL_MAX, point_dist = DBL_MAX, global_cent_dist, global_point_dist, temp_dist, rad_comp_dist;
	double *cent_comp_arr;

	cent_comp_arr = (double *)malloc(sizeof(double) * k);

	for (j = 0; j < q; j++) { // for each query
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
        //printf("%d has min cent dist %f idx at %d\n", world_rank, cent_dist, min_clust_idx);

		//calculate distance to each point in the closest cluster and get the local minimum
        loopstart = cluster_start[min_clust_idx];
        loopend = cluster_start[min_clust_idx]+cluster_size[min_clust_idx]*dim;

		for (l = loopstart; l < loopend; l+=dim) {
			//temp_dist = distance(dim, j, query, l, data);
            temp_dist = pnt2pntDistance(dim, j, query, l, data);
			count++;
			point_dist = MIN(point_dist, temp_dist);
			if (point_dist == temp_dist) {
				min_point_idx = l;
			}	
		}

        //printf("%d had min point_dist %f\n", world_rank, point_dist);

        //all procs need to have the same minimum distance to a point in the cluster
        MPI_Allreduce(&point_dist, &global_point_dist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        point_dist = global_point_dist;

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

        //all procs need to have the same minimum distance to a point in the cluster
        MPI_Allreduce(&point_dist, &global_point_dist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        //if we had the shortest distance, we will add our result to the result array
        if (point_dist == global_point_dist) {
		    for (i = 0; i < dim; i++) {
		    	result_pt[j*dim+i] = data[min_point_idx+i];
	    	}
        }
        
        point_dist = DBL_MAX;
        cent_dist = DBL_MAX;
	}

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

    if (world_rank == 0) 
        printf("Adjusting centroids\n...\n");

    //do the process until none of the data points change cluster (has a hardcoded limit of 50, can be modified or removed)
	while ((!doneflag) && (stallflag < 50)) {

		MPI_Barrier(MPI_COMM_WORLD);

		if (stallflag % 10 == 0) { //check for back and forthness
			for (i = 0; i < ndata; i++) {
				cycle_check_assign[i] = cluster_assign[i];
			}
		}

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
		stallflag++;
		cycle_cnt++;
	}

    if (world_rank == 0) { //printf global cluster sizes
        printf("\nCluster sizes: ");
        for (i = 0; i < k; i++) {
            printf("%d ", temp_clust_size[i]);
        }
        printf("\n");
    }

	for (i = 0; i < k; i++) {
		loopstart = cluster_start[i]*dim;
		if (i == k - 1) {
			loopend =  ndata*dim;
		}
		else {
			loopend = cluster_start[i+1]*dim;
		}
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
		loopstart = cluster_start[i]*dim;
		if (i == k - 1) {
			loopend = ndata*dim;
		}
		else {
			loopend = cluster_start[i+1]*dim;
		}
		//for each centroid, getting the sum of all points, each of their dimensions, in temp_centroid
		for (j = loopstart; j < loopend; j+=dim) {
			for (l = 0; l < dim; l++) {	
				temp_centroid[l] += data[j+l];
			}
		} 

		MPI_Allreduce(&temp_centroid[0], &temp_proc_cent[0], dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);

		//now get the averages for each dimension in temp_centroid and put it in cluster_centroid
		//printf("cluster centroid ");
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

/*double calculateClusterDistance(int dim, int data_idx, double *data, int clust_idx, double **cluster_centroid) {
	int i;
	double distance = 0.0;

	for (i = 0; i < dim; i++) {
		distance += pow(cluster_centroid[clust_idx][i] - data[data_idx+i], 2);
	}

	distance = sqrt(distance);

	return distance;
}*/

/*double calculatePointDistance(int dim, int clust_idx, int query_idx, double *query, double **cluster_centroid) {
	int i;
	double distance = 0.0;

	for (i = 0; i < dim; i++) {
		distance += pow(cluster_centroid[clust_idx][i] - query[query_idx+i], 2);
	}

    distance = sqrt(distance);

	return distance;
}

double distance(int dim, int query_idx, double *query, int data_idx, double *data) {
	int i;
	double distance = 0.0;

	for (i = 0; i < dim; i++) {
		distance += pow(data[data_idx+i] - query[query_idx+i], 2);
	}

    distance = sqrt(distance);

	return distance;
}*/

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
			//temp_distance = calculateClusterDistance(dim, i, data, j, cluster_centroid);
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
	//printf("cluster size, ");
	for (l = 0; l < k; l++) {
		cluster_start[l] = clust_size_idx_cnt;
		clust_size_idx_cnt += cluster_size[l];
		//printf("%d,", cluster_size[l]); 
	}
	//printf("\n");

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
	//printf("here1\n");

	for (j = 0; j < ndata*dim; j+=dim) { //for each data point
		for (i = 0; i < m; i++) { //for given m centroids
			if (i != place_idx) {
                distance = pnt2centDistance(dim, i, j, data, cluster_centroid);
				/*for (l = 0; l < dim; l++) {
					distance += pow(cluster_centroid[i][l] - data[j+l], 2); 		
				}	
				distance = sqrt(distance);*/
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
}




int main(int argc, char** argv) {
	

	MPI_Init(NULL, NULL);

	int kcheck = 1, rint, cycles, pointcnt, strsize, q;
    int k, dim, ndata, data_remainder;
    char *path;
	double r;
	double *proc_data, **cluster_centroid, *cluster_radius, *query, *result;
	int *cluster_assign, *cluster_size,  *cluster_start;
    
    strsize = strlen(argv[1]);
    k = atoi(argv[2]);
    dim = atoi(argv[3]);
    ndata = atoi(argv[4]);
    q = atoi(argv[5]);

    //get path from command line
    path = (char *)malloc(sizeof(char) * strsize);
    path[0] = '\0';
    strcat(path, argv[1]);
    
    /*FILE *file = fopen(path, "rb");

    if (file == NULL) { //if no file
        printf("File not found, program exited with code 0\n");
        return 0;
    }*/

    //seed random number
	srand(23);
	
	cluster_size = (int *)malloc(sizeof(int) * k);
	cluster_start = (int *)malloc(sizeof(int) * k);
	cluster_centroid = (double **)malloc(sizeof(double *) * k);
	cluster_radius = (double *)malloc(sizeof(double) * k);

	query = (double *)malloc(sizeof(double) * q * dim);
	result = (double *)malloc(sizeof(double) * q * dim);

	int i, j, world_size, world_rank, proc_chunk_size;
	float *ft_data;
    unsigned char *inc_data;

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
	inc_data = (unsigned char *)malloc(sizeof(unsigned char) * proc_chunk_size * dim);


	//initialize 2d arrays
	for (i = 0; i < k; i++) {
		cluster_centroid[i] = (double *)malloc(sizeof(double) * dim);
	}

	for (i = 0; i < k; i++) {
		cluster_size[i] = 0;
	}

    for (i = 0; i < q * dim; i++) {
        query[i] = 0.0;
    }
	
    //start reading binary file
    /*fseek(file, world_rank * proc_chunk_size * dim * sizeof(float), SEEK_SET);

    fread(ft_data, sizeof(float), proc_chunk_size*dim, file);*/
    readFloatBin(path, ft_data, proc_chunk_size * dim, world_rank);
    //convert float values into doubles because computation numbers can get quite large
    for (i = 0; i < proc_chunk_size * dim; i++) {
        proc_data[i] = (double) ft_data[i];
    }

    //initialize random query points unless they are coming from a file
	for (i = 0; i < q * dim; i+=1000) {
        rint = rand() % 1000; //for size purposes, only have 1 out of every 1000 genes be expressed only up to an amount of 20
		r = ((double) rand() / (double)(RAND_MAX)) * 20.0;
        
		query[rint] = r; 
		//printf("%f,", r);
	}

	printf("\n"); //the first 1
	for (i = 0; i < dim; i+=dim) {
		printf("%d) ", proc_chunk_size*world_rank+i/dim);
		for (j = 0; j < 10; j++) {
			printf("%f,", proc_data[i+j]);
		}
		printf("\n\n");
	}

	if (world_rank == 0) {
        //initialize random query points unless they are coming from a file
	    for (i = 0; i < q * dim; i+=1000) {
            rint = rand() % 1000; //for size purposes, only have 1 out of every 1000 genes be expressed only up to an amount of 20
	    	r = ((double) rand() / (double)(RAND_MAX)) * 20.0;
        
	    	query[rint] = r; 
		    printf("%f,", r);
    	}
        printf("\n");
        MPI_Bcast(&query[0], dim * q, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        /*printf("Query: ");
        for (i = 0; i < 1000; i++) {
            printf("%f, ", query[i]);
        }
        printf("\n");*/

		//initialize centroids
		rint = rand() % proc_chunk_size;
		for (i = 0; i < dim; i++) {
			cluster_centroid[0][i] = proc_data[rint*dim+i];
		}
	
		MPI_Bcast(&cluster_centroid[0][0], dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        printf("Now initializing centroids\n...\n");
	}

	else { // not root
        MPI_Bcast(&query[0], dim * q, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&cluster_centroid[0][0], dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

    //before while loops good to make sure everyone is on the same page
	MPI_Barrier(MPI_COMM_WORLD);


	while (kcheck < k) {
        //printf("%d kcheck %d\n", world_rank, kcheck);
		initializeCentroids(dim, proc_chunk_size, proc_data, k, cluster_centroid, kcheck, kcheck, world_rank, world_size);
		kcheck++;
		MPI_Barrier(MPI_COMM_WORLD);
	}

    if (world_rank == 0) {
    	/*printf("\nCentroids\n");
    	for (i = 0; i < k; i++) {
    		for (j = 0; j < 10; j++) {
	    		printf("%f,", cluster_centroid[i][j]);
	    	}
	    	printf("\n");
	    }*/
    }

	cycles = kmeans(dim, proc_chunk_size, proc_data, k, cluster_size, cluster_start, cluster_radius, cluster_centroid, cluster_assign, world_size, world_rank);
	
    /*printf("\nCentroid sizes ");
	for (i = 0; i < k; i++) {
		printf("%d,", cluster_size[i]);
	}
	printf("\n");*/


if (world_rank == 0) {

	printf("\ncycles %d\n", cycles);

    /*printf("Cluster radii: ");
    for (i = 0; i < k; i++) {
        printf("%f, ", cluster_radius[i]);
    }
    printf("\n");*/

	//printf("\n");
}
	pointcnt = search_kmeans(dim, ndata, proc_data, k, cluster_size, cluster_start, cluster_radius, cluster_centroid, query, result, world_size, world_rank, q);

	printf("number of points checked %d\n\n", pointcnt/q);

    free(ft_data);
	free(result);
	free(query);
	free(cluster_radius);
	free(cluster_centroid);
	free(cluster_start);
	free(cluster_size);

    free(proc_data);
    free(inc_data);

	MPI_Finalize();

	return 0;
}

