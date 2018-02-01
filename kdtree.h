/*
 *@param dim: dimension of dataset
 *@param ndata: the number of data points
 *@param data: array of the data points. Size ndata * dim
 *@param k: number of clusters that will be made
 *@param cluster_size: array of the sizes of each cluster. Size k
 *@param cluster_start: array of the starting points of each cluster in the data array. Size k
 *@param cluster_bndry: array of the boundries of each cluster. Size (2 * dim) * k
 *@param cluster_assign: array of each data point's cluster assignment. Size ndata
 *@param query: array holding the query points. Size number_of_queries(called Q) * dim
 *@param result_pt: array holding the closest point to each corresponding query point. Size Q * dim
 */

void bipartition(int dim, int i0, int im, double *data, int *cluster_size, int *cluster_start, double *cluster_bdry, double *cluster_centroid, int *cluster_assign);

void biparttracker(int dim, int ndata, int depth, int cluster, int i0, int im, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_bdry, double **cluster_centroid, int *cluster_assign);

void kdtree(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_bdry, double **cluster_centroid, int *cluster_assign);

int search_kdtree(int dim, int ndata, double *data, int k, int q, int *cluster_size, int *cluster_start, double **cluster_bdry, double *query, double *result_pt);

void runKDTree(char *path, int ndata, int dim, int k, int q, double *query);
