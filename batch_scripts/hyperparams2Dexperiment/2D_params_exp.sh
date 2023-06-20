: '
ARGUMENTS:
    - $1: path to the directory containing the files
    - $2: path to the directory where the segway runs will be stored
    - $3: path to the directory where the reproducibility results will be stored
    - $4: segment transition weighting START
    - $5: segment transition weighting END
    - $6: segment transition weighting STEP
    - $7: resolution (bin size) START
    - $8: resolution (bin size) END
    - $9: resolution (bin size) STEP
'

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
for seg_trans_wt in {$4..$5..$6}; do
    for res in {$7..$8..$9}; do
        DIRS_SUFFIX="seg_trans_wt-${seg_trans_wt}_res-${res}"
        python run.py /home/wxa38/scratch/files /home/wxa38/scratch/segway_runs/${DIRS_SUFFIX} /home/wxa38/scratch/reprod_results${DIRS_SUFFIX}