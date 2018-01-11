#define LABEL_DIM 33

void txtToFloatBin(char *read_path, char *write_path, float *data, int features, int ndata, int labels_flag);
void txtToCharBin(char *read_path, char *write_path, unsigned char *data, int features, int ndata, int labels_flag);
void readFloatBin(char *path, float *data, int chunk_size, int rank);
void readCharBin(char *path, unsigned char *data, int chunk_size, int rank);
