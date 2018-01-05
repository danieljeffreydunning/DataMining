#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dataFunctions.h"

#define LABEL_DIM 33

void txtToFloatBin(char *read_path, char *write_path, float *data, int features, int ndata, int labels_flag) {
	FILE* txtfile = fopen(read_path, "r");
    FILE* binfile = fopen(write_path, "wb");
 
	int i, j, k;
    char labs[LABEL_DIM];

    if (labels_flag == 1) {
        /*First we read the text file*/
        k = 0;
    	for (i = 0; i < ndata; i++) {
    		// Read the label
    	    fscanf(txtfile, "%f", &data[k]);
    		k++;
        
	    	// Read the features
        	for (j = 0; j < features; j++) {
    	    	fscanf(txtfile, "%f", &data[k]);
    	   	    k++;
    	    }
        }

        /* Now write to binary file */
        fwrite(&data[0], features * ndata, sizeof(float), binfile);
    }
    else {
        /*First we read the text file*/
        k = 0;
    	for (i = 0; i < ndata; i++) {
    		// Read the label
    	    fscanf(txtfile, "%s", labs);
        
	    	// Read the features
        	for (j = 0; j < features; j++) {
    	    	fscanf(txtfile, "%f", &data[k]);
    	   	    k++;
    	    }
        }

        /* Now write to binary file */
        fwrite(&data[0], features * ndata, sizeof(float), binfile);
    }

	fclose(binfile);
    fclose(txtfile);
}

void txtToCharBin(char *read_path, char *write_path, unsigned char *data, int features, int ndata, int labels_flag) {
	FILE* txtfile = fopen(read_path, "r");
    FILE* binfile = fopen(write_path, "wb");
 
	int i, j, k;
    char labs[LABEL_DIM];

    if (labels_flag == 1) {
        /*First we read the text file*/
        k = 0;
    	for (i = 0; i < ndata; i++) {
    		// Read the label
    	    fscanf(txtfile, "%s", &data[k]);
    		k++;
        
	    	// Read the features
        	for (j = 0; j < features; j++) {
    	    	fscanf(txtfile, "%s", &data[k]);
    	   	    k++;
    	    }
        }

        /* Now write to binary file */
        fwrite(&data[0], features * ndata, sizeof(unsigned char), binfile);
    }
    else {
        /*First we read the text file*/
        k = 0;
    	for (i = 0; i < ndata; i++) {
    		// Read the label
    	    fscanf(txtfile, "%s", labs);
        
	    	// Read the features
        	for (j = 0; j < features; j++) {
    	    	fscanf(txtfile, "%s", &data[k]);
    	   	    k++;
    	    }
        }

        /* Now write to binary file */
        fwrite(&data[0], features * ndata, sizeof(unsigned char), binfile);
    }

	fclose(binfile);
    fclose(txtfile);

}

void readFloatBin(char *path, float *data, int chunk_size) {
    FILE *file = fopen(path, "rb");

    if (file == NULL) { //if no file
        printf("File not found, program exited with code 0\n");
        return;
    }

    //starting position
    fseek(file, chunk_size * sizeof(float), SEEK_SET);

    //read from starting position with chunk size of chunk_size
    fread(data, sizeof(float), chunk_size, file);

    fclose(file);
}

void readCharBin(char *path, unsigned char *data, int chunk_size) {
    FILE *file = fopen(path, "rb");

    if (file == NULL) { //if no file
        printf("File not found, program exited with code 0\n");
        return;
    }

    //starting position
    fseek(file, chunk_size * sizeof(unsigned char), SEEK_SET);

    //read from starting position with chunk size of chunk_size
    fread(data, sizeof(unsigned char), chunk_size, file);

    fclose(file);
}
