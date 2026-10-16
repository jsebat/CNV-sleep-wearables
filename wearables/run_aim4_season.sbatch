#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --job-name=aim4_season
# Aim 4 models with season covariate (21, SEASON=1), all ages and 7-25.
export OPENBLAS_NUM_THREADS=1
A=/expanse/projects/sebat1/g2mh_data/empatica_1309_analysis; RS=$HOME/micromamba/envs/actig/bin/Rscript
echo "=== Aim 4 season, all ages"; SEASON=1 OUT_TAG=_season $RS $A/code/21_aim4_models.R $A
echo "=== Aim 4 season, 7-25"; SEASON=1 AGE_MIN=7 AGE_MAX=25 OUT_TAG=_season_7to25 $RS $A/code/21_aim4_models.R $A
