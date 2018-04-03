mpiexec -n 1 ./a.out ../data/ethylene_CO.bin -a 3 -nd 4208262 -d 16 -q 1000
mpiexec -n 1 ./a.out ../data/match_t_p.bin -a 3 -nd 721 -d 60483 -q 1000
mpiexec -n 1 ./a.out ../data/ethylene_CO.bin -a 3 -nd 4208262 -d 16 -q 10000
mpiexec -n 1 ./a.out ../data/match_t_p.bin -a 3 -nd 721 -d 60483 -q 10000
