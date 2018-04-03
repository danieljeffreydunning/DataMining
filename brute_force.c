#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include "brute_force.h"
#include "util/dataFunctions.h"
#include "util/compFunctions.h"

void runBruteForce(char *path, int ndata, int dim, int q, double *query, double *result) {
    int i, j;
    int *closestIdx;
    float *ft_data;
    double minDistance = DBL_MAX, tempDistance;
    double *data;

    closestIdx = (int *)malloc(sizeof(int) * q);
    ft_data = (float *)malloc(sizeof(float) * ndata * dim);
    data = (double *)malloc(sizeof(double) * ndata * dim);

    //start reading binary file
    readFloatBin(path, ft_data, ndata * dim, 0);
    //convert float values into doubles because computation numbers can get quite large
    for (i = 0; i < ndata * dim; i++) {
        data[i] = (double) ft_data[i];
    }

    for (i = 0; i < q*dim; i+=dim) {
        for (j = 0; j < ndata * dim; j+=dim) {
            tempDistance = pnt2pntDistance(dim, i, query, j, data);
            minDistance = MIN(minDistance, tempDistance);
            if (minDistance == tempDistance) {
                closestIdx[i/dim] = j;
            }
        }
        result[i/dim] = minDistance;
        minDistance = DBL_MAX;
    }

    for (i = 0; i < q; i++) {
        printf("%f\n", result[i]);
    }
    //printf("\n");

    free(closestIdx);
    free(ft_data);
    free(data);
    free(result);
}
