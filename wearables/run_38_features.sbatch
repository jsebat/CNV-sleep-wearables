#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --job-name=aim2_feat
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
A=/expanse/projects/sebat1/g2mh_data/empatica_1309_analysis
$HOME/micromamba/envs/actig/bin/python -I $A/code/38_aim2_features.py $A $A/aim2/features $SLURM_CPUS_PER_TASK
