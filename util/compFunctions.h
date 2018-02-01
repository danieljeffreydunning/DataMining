#ifndef COMPFUNCTIONS_H
#define COMPFUNCTIONS_H

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define LOG2(X) log((X)) / log(2)

/*
 * 
 * */

double pnt2pntDistance(int dim, int query_idx, double *query, int data_idx, double *data);
double pnt2centDistance(int dim, int clust_idx, int query_idx, double *query, double **cluster_centroid);
int dotprod(int dim, double *data, double **r, int d_idx, int r_idx);
double pnt2bdry(int dim, double **cluster_bdry, double *query, int qidx, int clust_idx);


#endif /*COMPFUNCTIONS_H*/
