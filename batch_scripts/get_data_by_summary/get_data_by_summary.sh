#!/bin/bash
#SBATCH --account=def-maxwl
#SBATCH --cpus-per-task=8
#SBATCH --mem=125G
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wxa38@sfu.ca
#SBATCH --output=get_data_by_summary_$(date).log

cd /home/wxa38/projects/def-maxwl/wxa38/SAGAconf
# conda activate SAGAconf
python get_data.py