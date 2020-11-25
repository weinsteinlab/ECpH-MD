#!/bin/sh
#SBATCH --job-name=bb6
#SBATCH -p dcs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --mem=40G
##SBATCH --mail-user=kots.katya@gmail.com
##SBATCH --mail-type=END

source ~/.bashrc

module load gcc/8.1.0/1
module load cuda/10.1

conda activate openmm_7.4.0 

python -u Exchange-min-replica.py
