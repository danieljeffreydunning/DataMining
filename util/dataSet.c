#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char** argv) {
    int i, j, seed, ndata, rintq, rintd, rintsign, q, dim = 2;
    double r;
    double *query;
    float *data;
    FILE* binfile = fopen("../../data/localDataSet.bin", "wb");

    seed = atoi(argv[1]);
    ndata = atoi(argv[2]);
    q = atoi(argv[3]);

    data = (float *)malloc(sizeof(float) * dim * ndata);
    query = (double *)malloc(sizeof(double) * dim * q);

    srand(1994);
    for (i = 0; i < ndata * dim; i++) {
        r = ((float) rand() / (float)(RAND_MAX)) * 100.0;
        data[i] = r;
    }

    srand(seed);
    for (i = 0; i < q * dim; i++) {
        r = ((double) rand() / (double)(RAND_MAX)) * 100.0;
		query[i] = r; 
	}

    //for 5 points for each query
    for (i = 0; i < 5 * q; i++) {
        rintq = rand() % q;
        rintd = rand() % ndata;
        for (j = 0; j < dim; j++) {
            rintsign = rand() % 2; //random sign
            r = ((double) rand() / (double)(RAND_MAX)) * 3.0;
            data[rintd+j] = fabs(query[rintq+j] + (r * pow(-1, rintsign))); //each dim will be within 3 units of original dim
        }
    }

    fwrite(&data[0], dim * ndata, sizeof(float), binfile);

    /*for (i = 0; i < dim; i++) {
        printf("query %f, data %f\n", query[rintq+i], data[rintd+i]);
    }*/

    fclose(binfile);

    return 0;
}
