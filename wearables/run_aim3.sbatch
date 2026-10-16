#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --job-name=aim3
# usage: sbatch [--export=ALL,NREP=2] run_aim3.sbatch <predict_out_subdir> <n_perm> [assoc:1/0]
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
A=/expanse/projects/sebat1/g2mh_data/empatica_1309_analysis; E=$HOME/micromamba/envs/actig/bin
[ "${3:-1}" = 1 ] && $E/Rscript $A/code/41_aim3_assoc.R $A
$E/python -I $A/code/42_aim3_predict.py $A $A/aim3/$1 $SLURM_CPUS_PER_TASK $2
NO_MEDS=1 $E/python -I $A/code/42_aim3_predict.py $A $A/aim3/${1}_nomeds $SLURM_CPUS_PER_TASK $2
