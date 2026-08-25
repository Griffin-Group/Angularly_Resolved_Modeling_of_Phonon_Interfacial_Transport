#!/bin/bash
#SBATCH --job-name=al
#SBATCH --partition=dirac1
#SBATCH --account=dirac1
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=24
##SBATCH --cpus-per-task=1
#SBATCH --time=40:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mhussien@lbl.gov
#SBATCH --output=Job-%j.out

module purge

module load intel-oneapi-compilers/2023.1.0
module load intel-oneapi-mpi/2021.10.0
module load intel-oneapi-mkl/2023.2.0

module list

echo "The job id is: $SLURM_JOB_ID"
echo "Starting VASP at `date`"

mpirun -np $SLURM_NTASKS  ~/bin/vasp643_std

echo "Finished VASP at `date`"
exit 0
