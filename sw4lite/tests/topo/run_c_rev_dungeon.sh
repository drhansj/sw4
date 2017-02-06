export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#quad, cache
echo "1 MPI, 64 cores, 1 HT"
export OMP_NUM_THREADS=64
mpirun -n 1 ./sw4lite_c gaussianHill-rev.in

echo "2 MPI, 32 cores, 1 HT"
export OMP_NUM_THREADS=32
mpirun -n 2 ./sw4lite_c gaussianHill-rev.in

echo "4 MPI, 16 cores, 1 HT"
export OMP_NUM_THREADS=16
mpirun -n 4 ./sw4lite_c gaussianHill-rev.in

echo "4 MPI, 16 cores, 2 HT"
export OMP_NUM_THREADS=32
mpirun -n 4 ./sw4lite_c gaussianHill-rev.in



