/*
 *@param dim: dimension of each point in the data set
 *@param ndata: number of data points
 *@param data: array of the data points. Size ndata * dim
 *@param m: number of hashing vectors that will be used. Must be played with to get best results. Start small
 *@param r: array of the hashing vectors. Size m * dim
 *@param b: array of random offset values. Size m
 *@param w: number of chunks the data is divided into. Must be played with to get best results. Start large
 *@param num_clusters: current number of clusters at a given time
 *@param cluster_size: array with the size of each cluster. Size num_clusters
 *@param cluster_start: array with the starting point of each cluster in the data array. Size num_clusters
 *@param H: array with each data point's hash value. Size ndata * m
 *@param hash_vals: array that hold all the known unique hash values. (initialized to) size m * ndata
 */ 
int LSH(int dim, int ndata, double *data, int m, double **r, double *b, double w, int num_clusters, int *cluster_size, int *cluster_start, int **H, int *hash_vals);

/*
 *@param temp_hash: array of single hash value for simplicity of comparison. Size m
 *@param running_cnt: running count of the number of unique hash values that have been found. 
 *@param query: array of query points. Size Q * dim
 *@param result: array of the closest point to each respective query. Size Q * dim
 */
void local_search(int dim, int ndata, int q0, double *data, int *cluster_size, int *cluster_start, int m, int *hash_vals, int *temp_hash, int running_cnt, double *query, double *result);

/**/
int check_hash(int **H, int *hash_vals, int *clust_cnt, int idx, int m, int running_cnt, int *hash_assign); 

/**
 *temp_has: array set to the same as hash_vals. Will be used to do the reordering of values and then copied back into hash_vals. Size m * ndata
 */
void rearrange_data(double *data, int *cluster_size, int *cluster_start, int *hash_assign, int *hash_vals, int **H, int running_cnt, int m, int ndata, int dim); 

void runLSH(char *path, int ndata, int dim, int m, int w, int q);
