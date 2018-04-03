#!/bin/bash
N=$1
A=$2
ND=$3
D=$4
Q=$5

FILE="../data/worldcitiespop.bin"

mpicc main.c kdtree.c kmeans.c bisecting_kmeans.c lsh.c brute_force.c util/dataFunctions.c util/compFunctions.c -lm

if [ "$A" -lt 2 ]; then
    K=$6
    mpiexec -n $N ./a.out $FILE -a $A -nd $ND -d $D -k $K -q $Q    
elif [ "$A" -eq 2 ]; then
    M=$6
    W=$7
    mpiexec -n $N ./a.out $FILE -a $A -nd $ND -d $D -m $M -w $W -q $Q
else
    mpiexec -n $N ./a.out $FILE -a $A -nd $ND -d $D -q $Q    
fi
