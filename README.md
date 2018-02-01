# DataMining
Projects for Master's project

All projects are run through main.c. 
To compile:
mpicc main.c kdtree.c kmeans.c lsh.c util/dataFunctions.c util/compFunctions.c -lm

## KDTree
To run:

mpiexec -n 1 ./a.out {*path to input bin file*} -a 0 -nd {*number of data points*} -d {*number of dimensions*} -k {*number of clusters*} -q {*number of query points*}  

## KMeans 
Simple serial version available. Main project utilizes MPI and reads data from a binary file. Currently only suited for float values in the binary file(which will be converted to doubles for computation. More options on the TODO list

To run:

mpiexec -n {*number of procs*} ./a.out {*path to input bin file*} -a 1 -nd {*number of data points*} -d {*number of dimensions*} -k {*number of clusters*} -q {*number of query points*}  


## LSH

To run:

mpiexec -n 1 ./a.out {*path to input bin file*} -a 2 -nd {*number of data points*} -d {*number of dimensions*} -m {*number hash values (start low)*} -w {*width of each segment (start high)*} -q {*number of query points*}  

## Bisecting KMeans

