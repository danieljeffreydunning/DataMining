#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dataFunctions.h"

int main(int argc, char** argv) {
    int strsize1, strsize2, dim, ndata;
    char *pathfrom, *pathto;
    float *data;

    

    strsize1 = strlen(argv[1]);
    strsize2 = strlen(argv[2]);
    dim = atoi(argv[3]);
    ndata = atoi(argv[4]);

    data = (float *)malloc(sizeof(float) * dim * ndata);
    pathfrom = (char *)malloc(sizeof(char) * strsize1);
    pathto = (char *)malloc(sizeof(char) * strsize2);

    pathfrom[0] = '\0';
    pathto[0] = '\0';

    strcat(pathfrom, argv[1]);
    strcat(pathto, argv[2]);

    //printf("%s\n", pathfrom);
    //printf("%s\n", pathto);
    txtToFloatBin(pathfrom, pathto, data, dim, ndata, 0); 

    return 0;
}
