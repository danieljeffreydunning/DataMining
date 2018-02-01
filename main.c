#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "util/dataFunctions.h"
#include "util/compFunctions.h"
#include "kdtree.h"
#include "kmeans.h"
#include "lsh.h"
//
int parseArgs(int argc, char** argv, int *alg, int *ndata, int *dim, int *k, int *m, int *w, int *q) {
    int i;
    char this_alg[] = "-a", this_ndata[] = "-nd", this_dim[] = "-d", this_k[] = "-k", this_m[] = "-m", this_w[] = "-w", this_q[] = "-q";
    
    for (i = 2; i < argc; i+=2) {
        if (strcmp(argv[i], this_alg) == 0) {
            *alg = atoi(argv[i+1]);    
        }
        else if (strcmp(argv[i], this_ndata) == 0) {
            *ndata = atoi(argv[i+1]);
        }
        else if (strcmp(argv[i], this_dim) == 0) {
            *dim = atoi(argv[i+1]);
        }
        else if (strcmp(argv[i], this_k) == 0) {
            *k = atoi(argv[i+1]);
        }
        else if (strcmp(argv[i], this_m) == 0) {
            *m = atoi(argv[i+1]);
        }
        else if (strcmp(argv[i], this_w) == 0) {
            *w = atoi(argv[i+1]);
        }
        else if (strcmp(argv[i], this_q) == 0) {
            *q = atoi(argv[i+1]);
        }
        else {
            printf("%s not equal to something\n", argv[i]);
            return 0;
        }
    }
    return 1;
}

int main(int argc, char** argv) {
    int alg, q, ndata, dim, k, m, w, strsize, nitems, i;
    char *filepath;
    double r;
    double *query;

    strsize = strlen(argv[1]);

    filepath = (char *)malloc(sizeof(char) * strsize);
    filepath[0] = '\0';
    strcat(filepath, argv[1]);

    //if no file found
    if (filepath == NULL) {
        printf("File not found, program will use default data set\n");
        return 0;
    }

    if (!parseArgs(argc, argv, &alg, &ndata, &dim, &k, &m, &w, &q)) {
        printf("An argument was entered that is not recognized. Goodbye\n");
        return 0;
    }

    query = (double *)malloc(sizeof(double) * dim * q); 

    //initialize random query points
    srand(23);
    for (i = 0; i < q * dim; i++) {
		r = ((double) rand() / (double)(RAND_MAX)) * 100.0;
		query[i] = r; 
	}

    //run the specified algorithm
    if (alg == 0) {
        runKDTree(filepath, ndata, dim, k, q, query);
    }
    else if (alg == 1) {
        MPI_Init(&argc, &argv);
        runKMeans(filepath, ndata, dim, k, q, query);
    }
    else if (alg == 2) {
        runLSH(filepath, ndata, dim, m, w, q, query);
    }
    else {
        printf("That is not a valid algorithm. goodbye\n");
    }

    free(filepath);
    free(query);

    return 0;
}
