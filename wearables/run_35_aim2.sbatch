#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --job-name=aim2_xgb
# usage: sbatch [--export=ALL,NREP=2] run_35_aim2.sbatch <out_subdir> <n_perm>
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
A=/expanse/projects/sebat1/g2mh_data/empatica_1309_analysis
$HOME/micromamba/envs/actig/bin/python -I $A/code/35_aim2_xgb.py $A/aim4/master/G2MH_wearable_phenotypes_per_participant.csv $A/aim2/$1 $SLURM_CPUS_PER_TASK $2
