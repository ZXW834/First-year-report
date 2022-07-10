#!/bin/sh


#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1  

echo Hello from Slurm
module load intel-compilers/2021.2.0
module load R
R CMD BATCH summarylogistic_with_informative_prior.R
