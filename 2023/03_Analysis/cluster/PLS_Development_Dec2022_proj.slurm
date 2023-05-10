#!/bin/bash
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH --gres=gpu:1
#SBATCH --output=/home/marta.vidalgarcia/Development/slurm-PLS_proliferation-shape_Dec2022_proj.out

# Load R
module load R

# Run your script
Rscript Rscript_cluster_PLS_proliferation_proj_corrected.R
