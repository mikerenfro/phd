#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --time=24:00:00

module load warp3d
module load miniconda
source activate my_anaconda
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${OMP_NUM_THREADS}
python plate_runner.py
