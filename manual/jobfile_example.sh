#!/bin/bash
#SBATCH --job-name=d2
#SBATCH --mail-type=ALL
#SBATCH --partition=global
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=alessandro.leonardi.ing@gmail.com
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=64M
#SBATCH --output=./results/d2/output.txt
#SBATCH --error=./results/d2/error.txt

export OMP_NUM_THREADS=4

./hybird.exe -c ./erosion2D_particles.cfg -d ./results -n d2
