/*
 *@param dim: dimension of dataset
 *@param ndata: the number of data points
 *@param data: array of the data points. Size ndata * dim
 *@param k: number of clusters that will be made
 *@param cluster_size: array of the sizes of each cluster. Size k
 *@param cluster_start: array of the starting points of each cluster in the data array. Size k
 *@param cluster_radius: array of radii for each cluster. Size k * dim
 *@param cluster_centroid. array of centroids for each cluster. Size k * dim
 *@param cluster_assign: array of each data point's cluster assignment. Size ndata
 *@param query: array holding the query points. Size number_of_queries(called Q) * dim
 *@param result_pt: array holding the closest point to each corresponding query point. Size Q * dim
*/

void calculateInitCentroids(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int clust1, int clust2);

int assignInitData(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_centroid, int *cluster_assign, int clust1, int clust2);

int bkInitCentroids(int dim, int ndata, double *data, int *cluster_size, int *cluster_start, double **cluster_centroid, int place_idx, int cent_idx, int cent_count, int *cluster_assign);

void runBKmeans(int dim, int k, int ndata, double *data, double **cluster_centroid);

