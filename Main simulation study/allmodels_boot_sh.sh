#!/bin/bash
#
#SBATCH --job-name=zdunias_allmodels_boot_batch_4_repeat_double
#SBATCH --output="./data/out/output_allmodels_%a_boot_batch_4_repeat_double.txt"
#
#SBATCH --ntasks=1
#SBATCH --time=07:30:00
#SBATCH --mem-per-cpu=20
#SBATCH --chdir="./"
#SBATCH --mail-user=z.s.dunias@uu.nl
#SBATCH --mail-type=ALL
#
#SBATCH --array=375

echo job is started in Thesis
srun Rscript "./Rcode/Execute_boot.R" $SLURM_ARRAY_TASK_ID 0 25 4 1
echo job is finished in Thesis
