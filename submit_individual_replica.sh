#!/bin/sh
#SBATCH --job-name=bb6
#SBATCH -p dcs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --mem=40G
##SBATCH --mail-user=CWID@med.cornell.edu
##SBATCH --mail-type=END

source ~/.bashrc

module load gcc/8.1.0/1
module load cuda/10.1

conda activate openmm_7.4.0


pH=$1
#iteration=$2
#nsteps=$3
#python -u run_replica.py $pH $iteration $nsteps
python -u run_replica.py $pH
