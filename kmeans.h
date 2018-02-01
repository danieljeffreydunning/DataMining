#include <mpi.h>

/*
 *@param dim: dimension of dataset
 *@param ndata: the number of data points
 *@param data: array of the data points. Size ndata * dim
 *@param k: number of clusters that will be made
 *@param cluster_size: array of the sizes of each cluster. Size k
 *@param cluster_start: array of the starting points of each cluster in the data array. Size k
 *@param cluster_radius: array of radii of each cluster. Size k * dim
 *@param cluster_centroid: array of centroids for each cluster. Size k * dim
 *@param cluster_assign: array of each data point's cluster assignment. Size ndata
 *@param world_size: size of MPI communication world
 *@param world_rank: individual processor rank
 *@param query: array holding the query points. Size number_of_queries(called q) * dim
 *@param result_pt: array holding the closest point to each corresponding query point. Size q * dim
 *@param q: number of query point 
 */

int kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank);

void assignData(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank);

void calculateCentroids(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int world_size, int world_rank);

int search_kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, double *query, double *result_pt, int world_size, int world_rank, int q);

void initializeCentroids(int dim, int ndata, double *data, int k, double **cluster_centroid, int place_idx, int m, int world_rank, int world_size);


void runKMeans(char *path, int ndata, int dim, int k, int q, double *query);
