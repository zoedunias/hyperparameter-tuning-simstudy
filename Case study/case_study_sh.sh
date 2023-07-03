#!/bin/bash
#
#SBATCH --job-name=zdunias_case_study
#SBATCH --output="./output_case_study.txt"
#
#SBATCH --ntasks=1
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=20
#SBATCH --chdir="./"
#SBATCH --mail-user=Z.S.Dunias-2@umcutrecht.nl
#SBATCH --mail-type=ALL
#
#SBATCH --array=1


echo job is started in case_study
srun Rscript "./Thesis analysis empirical data V9 - PhD.R" 
echo job is finished in case_study
