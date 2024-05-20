#!/bin/bash
#SBATCH --job-name=fit_models
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=cpuq
#SBATCH --mem=100GB
#SBATCH --output=fit_models_%j.log
#SBATCH --error=fit_models_%j.stderr
#SBATCH --mail-user=m.logan@aims.gov.au
#SBATCH --mail-type=ALL
Rscript 20_fit_models.R
