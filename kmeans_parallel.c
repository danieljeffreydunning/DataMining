#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define N 32
#define K 4 //at least 2 
#define DIM 2
#define Q 1

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define LOG2(X) log((X)) / log(2)

//function prototypes
int kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank);
void calculateCentroids(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank);
void assignData(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank);
int search_kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, double *query, double *result_pt, int world_size, int world_rank);
double calculateClusterDistance(int dim, int data_idx, double *data, int clust_idx, double **cluster_centroid);
double calculatePointDistance(int dim, int clust_idx, int query_idx, double *query, double **cluster_centroid);
double distance(int dim, int query_idx, double *query, int data_idx, double *data);
void initializeCentroids(int dim, int ndata, double *data, int k, double **cluster_centroid, int start_idx, int m, int world_rank, int world_size);
//end prototypes

int search_kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, double *query, double *result_pt, int world_size, int world_rank) {
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


int kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank) {
	int i, j, doneflag = 0, cycle_cnt = 1, loopstart, loopend, stallflag = 0, empty_idx = k+5, cycle_check = 0, f = 0, diff = 0, proc_diff = 0;
	int *temp_assign, *cycle_check_assign, max_dist = 0.0, temp_max;
	double *temp_radius;	

	temp_assign = (int *)malloc(sizeof(int) * ndata);
	cycle_check_assign = (int *)malloc(sizeof(int) * ndata);
	temp_radius = (double *)malloc(sizeof(double) * k);

	//first call of assignData
	assignData(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign, world_size, world_rank);

	//see if there are any empty clusters
	for (i = 0; i < k; i++) {
		if (cluster_size[i] == 0) {
			empty_idx = i;
			break;
		}	
	}
	//printf("proc %d empty_idx %d\n", world_rank, empty_idx);

	//make sure there are no empty clusters
	while ((empty_idx < k) && (f < k)) {
		f++;
		initializeCentroids(dim, ndata, data, k, cluster_centroid, empty_idx, k, world_rank, world_size); 
		assignData(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign, world_size, world_rank);
		empty_idx = k+5;

		//see if there are any empty clusters
		for (i = 0; i < k; i++) {
			if (cluster_size[i] == 0) {
				empty_idx = i;
				break;
			}	
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}	
	f = 0;

	calculateCentroids(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign, world_size, world_rank);

	while ((!doneflag) && (stallflag < 50)) {
		MPI_Barrier(MPI_COMM_WORLD);
		//printf("proc %d while loop iter %d\n", world_rank, stallflag);

		if (stallflag % 10 == 0) { //check for back and forthness
			for (i = 0; i < ndata; i++) {
				cycle_check_assign[i] = cluster_assign[i];
			}
		}

	/*	printf("\n");
		for (i = 0; i < ndata; i++) {
			printf("%d,", cluster_assign[i]);
		}
		printf("\n");
		printf("\n");*/
		//if (world_rank == 0) {
		/*for (i = 0; i < ndata*dim; i+=dim) {
			printf("%d) ", ndata*world_rank+i/dim);
			for (j = 0; j < dim; j++) {
				printf("%f,", data[i+j]);
			}
			printf("\t%d\n", cluster_assign[i/dim]);
		}*/
		//}

		//initialize array for comparison
		for (i = 0; i < ndata; i++) {
			temp_assign[i] = cluster_assign[i];
		}
		//do "next" assignment. if nothing changes, we will be stopping
		assignData(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign, world_size, world_rank);
		empty_idx = k+5;

		//see if there are any empty clusters
		for (i = 0; i < k; i++) {
			if (cluster_size[i] == 0) {
				empty_idx = i;
				break;
			}	
		}
		//printf("empty_idx %d\n", empty_idx);

		while ((empty_idx < k) && (f < k)) {
			f++;
			initializeCentroids(dim, ndata, data, k, cluster_centroid, empty_idx, k, world_rank, world_size); 
			assignData(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign, world_size, world_rank);
			//printf("%d\n", cluster_size[world_rank]);
			empty_idx = k+5;

			//see if there are any empty clusters
			for (i = 0; i < k; i++) {
				if (cluster_size[i] == 0) {
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
		//printf ("PROC %d PROC_DIFF %d\n", world_rank, proc_diff);

		if (proc_diff > 0) {
			doneflag = 0;
			calculateCentroids(dim, ndata, data, k, cluster_size, cluster_start, cluster_centroid, cluster_assign, world_size, world_rank);
		//	printf("proc %d got here\n", world_rank);
		} 

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
		/*if (world_rank == 0) {
		printf("centroids\n");
		for (i = 0; i < k; i++) {
			for (j = 0; j < dim; j++) {
				printf("%f,", cluster_centroid[i][j]);
			}
			printf("\n");
		}
		}*/
		diff = 0;
		proc_diff = 0;
		stallflag++;
		cycle_cnt++;
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
			temp_max = calculateClusterDistance(dim, j, data, i, cluster_centroid);
			max_dist = MAX(max_dist, temp_max);
		}
		cluster_radius[i] = max_dist;
		
	}
	//reduce MAX operation for cluster radiuses
	MPI_Allreduce(&cluster_radius[0], &temp_radius[0], k, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	for (i = 0; i < k; i++) {
		cluster_radius[i] = temp_radius[i];
	}

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
	free(temp_assign);	
	free(cycle_check_assign);
	free(temp_radius);

	return cycle_cnt;
}

void calculateCentroids(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank) {
	//keep in mind cluster_start values are in terms of cluster_assign
	int i, j, l, loopstart, loopend;
	double temp_centroid[DIM] = {0.0}, temp_proc_cent[DIM] = {0.0}; 
	int *temp_clust_size; //actually the global cluster size

	temp_clust_size = (int *)malloc(sizeof(int) * k); 

	//initialize temp cluster size, used to get all of total cluster size
	//for (i = 0; i < k; i++) {
	//	temp_clust_size[i] = cluster_size[i];
	//}
	
	for (i = 0; i < k; i++) {
		loopstart = cluster_start[i]*dim;
		if (i == k - 1) {
			loopend = ndata*dim;
		}
		else {
			loopend = cluster_start[i+1]*dim;
		}
		//printf("loop start %d loop end %d \n", loopstart, loopend);
		//for each centroid, getting the sum of all points, each of their dimensions, in temp_centroid
		for (j = loopstart; j < loopend; j+=dim) {
			for (l = 0; l < dim; l++) {	
				temp_centroid[l] += data[j+l];
			}
		} 

		MPI_Allreduce(&temp_centroid[0], &temp_proc_cent[0], dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		/*if (world_rank == 0) {
			printf("total new centroid %d ", i); 
			for (j = 0; j < dim; j++) {
				printf("%f,", temp_proc_cent[j]); 
			}	
			printf("\n");

			printf("cluster sizes ");
			for (j = 0; j < k; j++) {
				printf("%d,", cluster_size[j]); 
			}
			printf("\n");
		}*/

		//MPI_Allreduce(&cluster_size[0], &temp_clust_size[0], k, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);

		//now get the averages for each dimension in temp_centroid and put it in cluster_centroid
		//printf("cluster centroid ");
		for (j = 0; j < dim; j++) {
			//printf("%f\t,", temp_centroid[j]); 
			cluster_centroid[i][j] = temp_proc_cent[j] / cluster_size[i];
			temp_centroid[j] = 0.0;
		}
		//printf("\n");
	}	

	free(temp_clust_size);
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

void assignData(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank) {
	int i, j, l, cent_idx, temp_data_idx = 0, clust_size_idx_cnt = 0, empty_idx = k+5;
	double temp_distance, distance = 10000.0;
	double *temp_data, *temp_assign;
	int *temp_clust_size;

	temp_data = (double *)malloc(sizeof(double) * ndata * dim);
	temp_assign = (double *)malloc(sizeof(double) * ndata);
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
			temp_distance = calculateClusterDistance(dim, i, data, j, cluster_centroid);

			distance = MIN(distance, temp_distance);

			if (distance == temp_distance)
				cent_idx = j;
		}

		cluster_assign[i/dim] = cent_idx;
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
	//	printf("%d,", cluster_size[l]); 
	}
	//printf("\n");

	//get the global cluster sizes
	MPI_Allreduce(&cluster_size[0], &temp_clust_size[0], k, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 

	for (i = 0; i < k; i++) {
		cluster_size[i] = temp_clust_size[i];
	}

	free(temp_clust_size);
	free(temp_assign);
	free(temp_data);
}

void initializeCentroids(int dim, int ndata, double *data, int k, double **cluster_centroid, int place_idx, int m, int world_rank, int world_size) {
	int i, j, l, centroid_point;
	double distance = 0.0, point_min = 10000.0, total_max = 0.0, global_max;
	//start_idx will be 0 in initialization case, but will be the cluster to be redone if calling this function for an empty cluster

	double *temp_point;

	temp_point = (double *)malloc(sizeof(double) * dim);
	

	for (j = 0; j < ndata*dim; j+=dim) { //for each data point
		for (i = 0; i < m; i++) { //for given m centroids
			if (i != place_idx) {
				for (l = 0; l < dim; l++) {
					distance += pow(cluster_centroid[i][l] - data[j+l], 2); 		
				}	
				distance = sqrt(distance);
				//get min distance to a centroid for that point
				point_min = MIN(point_min, distance);
				//printf("min distance for point %d is %f\n", j/dim, point_min);
			 }
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

	//save our temp point. temp point will be used after broadcast to store value of global centroid
	for (i = 0; i < dim; i++) {
		temp_point[i] = data[centroid_point+i];
	}

	MPI_Allreduce(&total_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	//if we were the owner of global max, we had the greatest and must broadcast our data point
	if (global_max == total_max) {
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

	int kcheck = 1, rint, cycles, pointcnt, s1, s2, s3, s4;
	double r;
	double *data, **cluster_centroid, *cluster_radius, *query, *result;
	int *cluster_assign, *cluster_size,  *cluster_start;

	cluster_assign = (int *)malloc(sizeof(int) * N);
	cluster_size = (int *)malloc(sizeof(int) * K);
	cluster_start = (int *)malloc(sizeof(int) * K);
	cluster_centroid = (double **)malloc(sizeof(double *) * K);
	cluster_radius = (double *)malloc(sizeof(double) * K);

	//query = (double *)malloc(sizeof(double) * Q * DIM);
	//result = (double *)malloc(sizeof(double) * DIM);

	int i, j, world_size, world_rank, proc_chunk_size;
	double *proc_data;

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	proc_chunk_size = N / world_size;

	proc_data = (double *)malloc(sizeof(double) * proc_chunk_size * DIM);


	//initialize 2d arrays
	for (i = 0; i < K; i++) {
		cluster_centroid[i] = (double *)malloc(sizeof(double) * DIM);
	}

	for (i = 0; i < K; i++) {
		cluster_size[i] = 0;
	}

	/*s1 = (int) argv[1];
	s2 = (int) argv[2];
	s3 = (int) argv[3];
	s4 = (int) argv[4];*/
/*	for (i = 0; i < N*DIM; i+=DIM) {
		td = calculateDistance(DIM, i, data, cluster_centroid); 
		d = MAX(d, td);
		if (td == d) {
			cent_idx = i;
		}
	}
	for (i = 0; i < dim; i++) {
		
	}*/
	
	

//if (world_rank == 0) {
	
	//data = (double *)malloc(sizeof(double) * N * DIM);

	srand(atoi(argv[world_rank+1]));
	
	for (i = 0; i < proc_chunk_size*DIM; i+=DIM) {
		for (j = 0; j < DIM; j++) {
			r = ((double)rand() / (double)(RAND_MAX)) * 100.0;
			proc_data[i+j] = r;
			//printf("%f,", r);
		}
		//printf("\n");
	}
	printf("\n");


	/*for (i = 0; i < Q * DIM; i++) {
		r = ((double)rand() / (double)(RAND_MAX)) * 100.0;
		query[i] = r; 
		printf("%f,", r);
	}*/


	printf("\n");
	for (i = 0; i < proc_chunk_size*DIM; i+=DIM) {
		printf("%d) ", proc_chunk_size*world_rank+i/DIM);
		for (j = 0; j < DIM; j++) {
			printf("%f,", proc_data[i+j]);
		}
		printf("\n");
	}
	printf("\n");

	/*printf("\nCentroid\n");
	for (i = 0; i < K; i++) {
		for (j = 0; j < DIM; j++) {
			printf("%f,", cluster_centroid[i][j]);
		}
		printf("\n");
	}*/
//	for (i = 0; i < proc_chunk_size*DIM; i++) {
//		proc_data[i] = data[i];
//	}

	/*printf("\nproc %d data\n", world_rank);
	for (i = 0; i < proc_chunk_size*DIM; i+=DIM) {
		printf("%d) ", i/DIM);
		for (j = 0; j < DIM; j++) {
			printf("%f,", proc_data[i+j]);
		}
		printf("\n");
	}*/

	/*for (i = 1; i < world_size; i++) {
		//send the portions of data
		MPI_Send(&data[i*proc_chunk_size*DIM], proc_chunk_size*DIM, MPI_DOUBLE, i, K+1, MPI_COMM_WORLD);	
		//send the centroids
	//	for (j = 0; j < K; j++) {
		//MPI_Send(&cluster_centroid[0][0], DIM, MPI_DOUBLE, i, 30, MPI_COMM_WORLD);
	//	}
	}*/

	if (world_rank == 0) {
		//initialize centroids
		rint = rand() % N;
		for (i = 0; i < DIM; i++) {
			cluster_centroid[0][i] = data[rint+i];
		}
	
		MPI_Bcast(&cluster_centroid[0][0], DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	//free(data);

//}// end of root proc initialization

	else { // not root
		//receive the chunk of data
		//MPI_Recv(&proc_data[0], proc_chunk_size*DIM, MPI_DOUBLE, 0, K+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		/*printf("\nproc %d data\n", world_rank);
		for (i = 0; i < proc_chunk_size*DIM; i+=DIM) {
			printf("%d) ", i/DIM);
			for (j = 0; j < DIM; j++) {
				printf("%f,", proc_data[i+j]);
			}
			printf("\n");
		}*/

		// receive each chunk of the centroid
		//for (i = 0; i < K; i++) {
		//MPI_Recv(&cluster_centroid[0][0], DIM, MPI_DOUBLE, 0, 30, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//}
		MPI_Bcast(&cluster_centroid[0][0], DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	}

	MPI_Barrier(MPI_COMM_WORLD);


	while (kcheck < K) {
		initializeCentroids(DIM, proc_chunk_size, proc_data, K, cluster_centroid, kcheck, kcheck, world_rank, world_size);
		kcheck++;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	/*if (world_rank < 2) {
		printf("\nproc %d centroids\n\n", world_rank);
		for (i = 0; i < K*DIM; i+=DIM) {
			printf("%d) ", i/DIM);
			for (j = 0; j < DIM; j++) {
				printf("%f,", cluster_centroid[i/DIM][j]);
			}
			printf("\n");
		}
		
	}*/
	//assignData(DIM, proc_chunk_size, proc_data, K, cluster_size, cluster_start, cluster_centroid, cluster_assign, world_size, world_rank);

	cycles = kmeans(DIM, proc_chunk_size, proc_data, K, cluster_size, cluster_start, cluster_radius, cluster_centroid, cluster_assign, world_size, world_rank);

if (world_rank == 0) {

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

	printf("\ncycles %d\n", cycles);

	printf("\n");
}
/*	pointcnt = search_kmeans(DIM, N, data, K, cluster_size, cluster_start, cluster_radius, cluster_centroid, query, result);

	for (i = 0; i < Q; i++) {
		for (j = 0; j < DIM; j++) {
			printf("%f,", result[i*DIM+j]);
		}
		printf("\n");
	}

	printf("number of points checked %d\n", pointcnt/Q);*/

	//free(result);
	//free(query);
	free(cluster_radius);
	free(cluster_centroid);
	free(cluster_start);
	free(cluster_size);

	MPI_Finalize();

	return 0;
}

