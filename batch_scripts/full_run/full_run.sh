#!/bin/bash
#SBATCH --account=def-maxwl
#SBATCH --cpus-per-task=24
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wxa38@sfu.ca
#SBATCH --output=full_run_00.log

cd /home/wxa38/projects/def-maxwl/wxa38/SAGAconf
# conda activate SAGAconf
python run.py /home/wxa38/scratch/files /home/wxa38/scratch/segway_runs /home/wxa38/scratch/reprod_results