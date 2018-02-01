#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "kdtree.h"
#include "util/dataFunctions.h"
#include "util/compFunctions.h"

void bipartition(int dim, int i0, int im, double *data, int *cluster_size, int *cluster_start, double *cluster_bdry, double *cluster_centroid, int *cluster_assign) {
    int loopstart = i0, loopend = im + 1, temp_cluster_size = im - i0 + 1;
    int *temp_assign;
    int i, j, clust1st = 0, clust2st = temp_cluster_size - 1, clust1sizecnt=0, clust2sizecnt=0;
    int tempclust1bound = clust1st, tempclust2bound = clust2st, clust1bnd = i0, clust2bnd = im; //in array, right bound for clust1 (left is start) and left bound for clust2 (end is right) both inclusive
    double *var_arr, *sum_arr, *mean_arr, *min_arr, *max_arr;
    double var_val = 0.0;
    int var_idx;
   
    temp_assign = (int *)malloc(sizeof(int) * temp_cluster_size);
    var_arr = (double *)malloc(sizeof(double) * dim);
    sum_arr = (double *)malloc(sizeof(double) * dim);
    mean_arr = (double *)malloc(sizeof(double) * dim);
    min_arr = (double *)malloc(sizeof(double) * dim);
    max_arr = (double *)malloc(sizeof(double) * dim);
    
    //initialize arr's to
    for (i = 0; i < dim; i++) {
        var_arr[i] = 0.0;
        sum_arr[i] = 0.0;
        mean_arr[i] = 0.0;
        max_arr[i] = 0.0;
        min_arr[i] = DBL_MAX;
    }
    
    for (i = loopstart, j = 0; i < loopend; i++, j++) {
        temp_assign[j] = cluster_assign[i];
    }
    
    //calculate dim sums
    for (i = 0; i < dim; i++) {
        for (j = loopstart; j < loopend; j++) { //j is index of cluster_assign
            sum_arr[i] += data[cluster_assign[j]+i];
            //also keep track of min and max values for boundries
            min_arr[i] = MIN(min_arr[i], data[cluster_assign[j]+i]);
            max_arr[i] = MAX(max_arr[i], data[cluster_assign[j]+i]);
        }
        //calculate dim means
        mean_arr[i] = sum_arr[i] / temp_cluster_size;
        
        //calculate dim variances
        for (j = loopstart; j < loopend; j++) {
            var_arr[i] += pow(data[cluster_assign[j]+i] - mean_arr[i], 2);
        }
    }
    
    //choose variance and save index
    for (i = 0; i < dim; i++) {
        var_val = MAX(var_val, var_arr[i]);
        //save index
        if (var_val == var_arr[i])
            var_idx = i;
    }

    //partition into each side of temp_assign based on the mean
    //cluster assign holds the values of the original data indexes for those points. used for rearranging at the end
    for (i = loopstart; i < loopend; i++) {
        //if the data point dim value at var_idx is less than the mean variance, it goes on the left side, else it goes on the right
        if (data[cluster_assign[i]+var_idx] < mean_arr[var_idx]) {
            temp_assign[clust1st] = cluster_assign[i];
            tempclust1bound = clust1st;
            clust1st++;
            clust1sizecnt++;
        }
        else {
            temp_assign[clust2st] = cluster_assign[i];
            tempclust2bound = clust2st;
            clust2st--;
            clust2sizecnt++;
        }
    }
    clust1bnd = clust1bnd + (clust1sizecnt - 1);
    clust2bnd = clust2bnd - (clust2sizecnt - 1);
    
    //reassign cluster_assign
    for (i = loopstart, j = 0; i < loopend; i++, j++) {
        cluster_assign[i] = temp_assign[j];
    }
    
    cluster_size[0] = clust1sizecnt;
    cluster_size[1] = clust2sizecnt;
    cluster_start[0] = i0;
    cluster_start[1] = clust2bnd;
    
    //update cluster boundries
    for (j = 0; j < 2; j++) { //each cluster
        for (i = 0; i < 2*dim; i+=2) { //each dim
            if (((i/2) == var_idx) && (j == 0)) { //1st cluster variance bounds
                cluster_bdry[j*2*dim+i] = min_arr[i/2];
                cluster_bdry[j*2*dim+i+1] = mean_arr[i/2];
                
            }
            else if (((i/2) == var_idx) && (j == 1)) { //2nd cluster variance bounds
                cluster_bdry[j*2*dim+i] = mean_arr[i/2];
                cluster_bdry[j*2*dim+i+1] = max_arr[i/2];
            }
            else { //either cluster general bounds
                cluster_bdry[j*2*dim+i] = min_arr[i/2];
                cluster_bdry[j*2*dim+i+1] = max_arr[i/2];
            }
        }
    }
    
    free(temp_assign);
    free(var_arr);
    free(sum_arr);
    free(mean_arr);
    free(min_arr);
    free(max_arr);
}

void biparttracker(int dim, int ndata, int depth, int cluster, int i0, int im, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_bdry, double **cluster_centroid, int *cluster_assign) {
    
    int i, j, clust1im, clust1i0, clust2i0, clust2im;
    int  *this_cluster_size,  *this_cluster_start;
    double *this_cluster_boundry, *this_cluster_centroid;
    
    if (depth > LOG2(k)) {
        //do nothing, reached max K
        return;
    }
    else if (depth == LOG2(k)) {
        /*for (i = cluster_start[cluster]; i < cluster_start[cluster] + cluster_size[cluster]; i++) {
            cluster_assign[i] = cluster;
        }*/
    }
    else {
        //allocate memory for single cluster
        this_cluster_size = (int *)malloc(sizeof(int) * 2);
        this_cluster_start = (int *)malloc(sizeof(int) * 2);
        this_cluster_boundry = (double *)malloc(sizeof(double) * 2 * 2 * dim);
        this_cluster_centroid = (double *)malloc(sizeof(double) * 2 * 2 * dim);
        
        
        //call biparttracker
        bipartition(dim, i0, im, data, this_cluster_size, this_cluster_start, this_cluster_boundry, this_cluster_centroid, cluster_assign);
        
        //simple calculation to keep a lot of arithmetic out of function call
        clust1im = this_cluster_start[0] + this_cluster_size[0] - 1;
        clust1i0 = i0;
        clust2im = this_cluster_start[1] + this_cluster_size[1] - 1;
        clust2i0 = this_cluster_start[1];
        
        //assign the global cluster starts/sizes
        for (i = 0; i < 2; i++) {
            if (depth == (LOG2(k) - 1)) {
                //cluster*2+i is the index of the next children
                printf("start %d size %d cluster %d\n", this_cluster_start[i], this_cluster_size[i], cluster*2+i);
                cluster_size[cluster*2+i] = this_cluster_size[i];
                cluster_start[cluster*2+i] = this_cluster_start[i];
            }
            for (j = 0; j < dim*2; j++) {
                cluster_bdry[cluster*2+i][j] = this_cluster_boundry[2*dim*i+j];
            }
        }

        //free allocated memory
        free(this_cluster_size);
        free(this_cluster_start);
        free(this_cluster_boundry);
        free(this_cluster_centroid);
        
        biparttracker(dim, ndata, depth+1, cluster*2, clust1i0, clust1im, data, k, cluster_size, cluster_start, cluster_bdry, cluster_centroid, cluster_assign);
        biparttracker(dim, ndata, depth+1, cluster*2+1, clust2i0, clust2im, data, k, cluster_size, cluster_start, cluster_bdry, cluster_centroid, cluster_assign);
        
    }
    
}

void kdtree(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_bdry, double **cluster_centroid, int *cluster_assign) {
    
    double *temp_data;
    int i, j, a, b;
    
    biparttracker(dim, ndata, 0, 0, 0, ndata-1, data, k, cluster_size, cluster_start, cluster_bdry, cluster_centroid, cluster_assign);
    //store data in temp data so we can rearrange data
    temp_data = (double *)malloc(sizeof(double) * ndata * dim);
    for (i = 0; i < ndata * dim; i++) {
        temp_data[i] = data[i];
    }
   
    //rearrange data 
    for (i = 0; i < ndata * dim; i+=dim) {
        a = cluster_assign[i/dim];
        for (j = 0; j < dim; j++) {
            data[i+j] = temp_data[a+j];
        }
    }

    //after kdtree is complete, need to make cluster start for data array, not cluster assign
    for (i = 0, j = 0; i < k; i++) {
        cluster_start[i] = j;
        j += cluster_size[i] * dim;
    }
    
    free(temp_data);
}

int search_kdtree(int dim, int ndata, double *data, int k, int q0, int *cluster_size, int *cluster_start, double **cluster_bdry, double *query, double *result_pt) {
    
    int i, j, q, qidx, min_clust_idx, point_idx, point_check_cnt = 0, loopstart, loopend;
    double min0, max0, calcdist, *clust_dist_arr, compare_dist = DBL_MAX, temp_dist = DBL_MAX;

    clust_dist_arr = (double *)malloc(sizeof(double) * k);
    
    for (q = 0; q < q0*dim; q+=dim) { //for each query
        qidx = q;
        for (i = 0; i < k; i++) { //for each cluster
            calcdist = pnt2bdry(dim, cluster_bdry, query, qidx, i);
            clust_dist_arr[i] = calcdist;
            
            compare_dist = MIN(compare_dist, calcdist); //find smallest distance to a cluster
            if (compare_dist == calcdist) {
                //we need the index of this cluster
                min_clust_idx = i;
            }
        }

        //keep "long" indexes out of for-loops
        loopstart = cluster_start[min_clust_idx];
        loopend = cluster_start[min_clust_idx] + cluster_size[min_clust_idx]*dim;
        
        //find closest point in the cluster we're in or closest to
        //compare dist now changes its meaning from distance to cluster to distance to points
        compare_dist = DBL_MAX;
        for (i = loopstart; i < loopend; i+=dim) {
            temp_dist = pnt2pntDistance(dim, qidx, query, i, data);
            compare_dist = MIN(temp_dist, compare_dist);
            if (compare_dist == temp_dist) {
                point_idx = i;
            }
            point_check_cnt++; //checked a new point
        }
        printf("closest point in cluster: %d with distance of %f\n", min_clust_idx, compare_dist);
        
        //check if any clusters are closer than the given point
        for (i = 0; i < k; i++) {
            if (i == min_clust_idx) {}//already checked this cluster
            else {
                loopstart = cluster_start[i];
                loopend = cluster_start[i] + cluster_size[i]*dim;
                if (clust_dist_arr[i] < compare_dist) { //cluster is closer than closest point in current cluster
                    for (j = loopstart; j < loopend; j+=dim) { //check all points in this cluster
                        temp_dist = pnt2pntDistance(dim, qidx, query, j, data);
                        compare_dist = MIN(temp_dist, compare_dist);
                        if (compare_dist == temp_dist) {
                            point_idx = i;
                            min_clust_idx = i;
                        }
                        point_check_cnt++;
                    }
                    
                }
            }
        }
        
        printf("[new] closest point in cluster: %d with a distance of %f\n", min_clust_idx, compare_dist);
        for (i = 0; i < dim; i++) {
            result_pt[qidx+i] = data[point_idx+i];
        }
        compare_dist = DBL_MAX;
        temp_dist = DBL_MAX;    
    }
    
    free(clust_dist_arr);
    
    return point_check_cnt / q0;
}

void runKDTree(char *path, int ndata, int dim, int k, int q, double *query) {
    int i, j, l, avg_checks, assign_idx;
    int *cluster_assign, *cluster_size,  *cluster_start;
    float *ft_data;
    double *data, **cluster_boundry, **cluster_centroid, *result;
    
    ft_data = (float *)malloc(sizeof(float) * ndata * dim);
    data = (double *)malloc(sizeof(double) * ndata * dim);
    cluster_assign = (int *)malloc(sizeof(int) * ndata);
    cluster_size = (int *)malloc(sizeof(int) * k);
    cluster_start = (int *)malloc(sizeof(int) * k);
    cluster_boundry = (double **)malloc(sizeof(double *) * k);
    cluster_centroid = (double **)malloc(sizeof(double *) * k);   
    result = (double *)malloc(sizeof(double) * q * dim);
    
    //initialize 2d arrays
    for (i = 0; i < k; i++) {
        cluster_boundry[i] = (double *)malloc(sizeof(double) * dim * 2);
        cluster_centroid[i] = (double *)malloc(sizeof(double) * dim);
    }
       
    //start reading binary file
    readFloatBin(path, ft_data, ndata * dim, 0);
    //convert float values into doubles because computation numbers can get quite large
    for (i = 0; i < ndata * dim; i++) {
        data[i] = (double) ft_data[i];
    }
    
    //initialize cluster_assign
    for (i = 0; i < ndata; i++) {
        cluster_assign[i] = i * dim;
    }
    
    kdtree(dim, ndata, data, k, cluster_size, cluster_start, cluster_boundry, cluster_centroid, cluster_assign);

    /*for (i = 0; i < k; i++) {
        printf("L%f R%f B%f T%f\n", cluster_boundry[i][0], cluster_boundry[i][1], cluster_boundry[i][2], cluster_boundry[i][3]);
    }*/
    for (i = 0; i < k; i++) {
        printf("%d %d\n", cluster_start[i], cluster_size[i]);
    }
    
    //give new values for cluster_assign for python visualization
    assign_idx = 0;
    for (i = 0; i < k; i++) {
        for (j = 0; j < cluster_size[i]; j++) {
            cluster_assign[assign_idx] = i;
            assign_idx++;
        }
    }

    write_results(dim, ndata, data, cluster_assign);

    /*for (i = 0; i < k; i++) {
        for (l = cluster_start[i]; l < cluster_start[i] + cluster_size[i]*dim; l+=dim) {
            printf("%d)", i);
            for (j = 0; j < dim; j++) {
                printf("\t%d\t%f\n", cluster_assign[l/dim], data[l+j]);
            }
        }
    }*/

    avg_checks = search_kdtree(dim, ndata, data, k, q, cluster_size, cluster_start, cluster_boundry, query, result);
    
    printf("Average number of checks was %d\n", avg_checks);
    printf("\n");
    
    for (i = 0; i < k; i++) {
        free(cluster_boundry[i]);
        free(cluster_centroid[i]);
    }
    
    free(data);
    free(cluster_size);
    free(cluster_start);
    free(cluster_boundry);
    free(cluster_centroid);
    free(cluster_assign);
    free(result);
}


