#!/bin/bash
#SBATCH --account=def-maxwl
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wxa38@sfu.ca
#SBATCH --output=rsync_CD14s.log

rsync -av ~/scratch/files/ ~/scratch/files