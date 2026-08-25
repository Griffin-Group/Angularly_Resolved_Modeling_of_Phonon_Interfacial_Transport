module load vasp/6.4.3-gpu

# One can use up to 16 OpenMP threads-per-MPI-rank when using
#  4 GPUs-per-node.
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Run 8 MPI ranks (-n 8) with 8 GPUs (-G 8) on 2 nodes:
#  The number of MPI ranks MUST match the number of GPUs.
#  The number of GPUs-per-node may not exceed 4.
#  Set -c ("--cpus-per-task") = 32 for the ideal MPI rank spacing.
srun -n 16 -c 32 -G 16 --cpu-bind=cores --gpu-bind=none vasp_std
