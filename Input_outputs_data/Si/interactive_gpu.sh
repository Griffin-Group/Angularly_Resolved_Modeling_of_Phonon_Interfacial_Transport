# Requesting an interactive session
#salloc --nodes=4 --qos=interactive --time=04:00:00 --constraint=gpu --gpus=16 --account=m3451_g

# Environment settings
export FI_CXI_DEFAULT_TX_SIZE=16384
export FI_CXI_RDZV_THRESHOLD=65536

# OpenMP Options
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# VASP Module
module load vasp/6.4.3-gpu
export PATH=$CFS/m3349/codes/vasp/vasp.6.4.3-gpu/bin/:$PATH

echo "The job id is: $SLURM_JOB_ID"
echo "Starting VASP at $(date)"

# Running VASP
srun --nodes=4 --gpus=16 --ntasks=16 --cpus-per-task=32 --ntasks-per-node=4 --gpus-per-task=1 --cpu-bind=cores --gpu-bind=none --unbuffered vasp_std &> output.log &

echo "Finished VASP at $(date)"
exit 0
