
module purge
module load intel/2023.2.0  
module load openmpi/5.0.7


# OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=true

export PATH=/pscratch/sd/m/mhussien/bin/EPW/q-e-EPW-6.0/bin:$PATH


srun -n 128 -c 2 --cpu-bind=cores q2r.x -in q2r.in > q2r.out
srun -n 128 -c 2 --cpu-bind=cores matdyn.x -in freq.in > freq.out
