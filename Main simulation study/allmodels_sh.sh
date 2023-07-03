#!/bin/bash
#
#SBATCH --job-name=zdunias_allmodels_batch_9
#SBATCH --output="./data/out/output_allmodels_%a_batch_9.txt"
#
#SBATCH --ntasks=1
#SBATCH --time=05:30:00
#SBATCH --mem-per-cpu=20
#SBATCH --chdir="./"
#SBATCH --mail-user=z.s.dunias@uu.nl
#SBATCH --mail-type=ALL
#
#SBATCH --array=1-275,301-450


echo job is started in Thesis
srun Rscript "./Rcode/Execute.R" $SLURM_ARRAY_TASK_ID 0 25 9 0
echo job is finished in Thesis
