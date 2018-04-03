#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    int comCnt = 0, i, j, dFlag = 0, ndata = 3173958;
    float var;
    float *data;
    char line[1000], *token, seps[] = ","; 
    static const char filename[] = "../../data/worldcitiespop.txt", binfilename[] = "../../data/worldcitiespop.bin";
    FILE *file = fopen (filename, "r");
    FILE *binfile = fopen(binfilename, "wb");

    data = (float *)malloc(sizeof(float) * ndata * 2);

    if (file != NULL) {
        while ( fgets ( line, sizeof line, file ) != NULL ) /* read a line */ {
            fputs (line, stdout); /* write the line */
            for (i = 0; i < 1000; i++) {
                if (line[i] == ',') {
                    comCnt++;
                    if (comCnt == 5) {
                        dFlag = 1;
                        break;
                    }
                }                    
            }

            if (dFlag == 1) {
                token = strtok(&line[i], seps);
                while (token != NULL) {
                    sscanf(token, "%f", &var);
                    data[j++] = var;
                    token = strtok(NULL, seps);
                }
            }

            dFlag = 0;
            comCnt = 0;
            memset(line,0,sizeof(line));
        }
        fclose (file);

        /*for (i = 0; i < 100; i++) {
            printf("%f, ", data[i]);
        }
        printf("\n\n");*/

    }
    fwrite(&data[0], ndata * 2, sizeof(float), binfile);

    free(data);

    return 0;
}
