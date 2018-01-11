# DataMining
Projects for Master's project

# KDTree

# KMeans 
Simple serial version available. Main project utilizes MPI and reads data from a binary file. Currently only suited for float values in the binary file(which will be converted to doubles for computation. More options on the TODO list

To comiple:
mpicc kmeans.c

To run
mpiexec -n {number of processors to be use}] ./a.out {path to input bin file} {number of clusters to use} {number of dimensions} {number of data points} {number of query points}     

# Bisecting KMeans

# LSH
