#!/bin/bash
#SBATCH -p q_intel_share
#SBATCH -n 4
#SBATCH -o job.out
 
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2

mpirun -n 4 ./direct_solver /online1/Solver25/solverchallenge25_01/solverchallenge25_01_A.mtx /online1/Solver25/solverchallenge25_01/solverchallenge25_01_b.rhs 1 0 0



