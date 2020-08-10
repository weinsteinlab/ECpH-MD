#!/bin/sh
#SBATCH --job-name=bb6
#SBATCH -p panda-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --mem=10G
##SBATCH --mail-user=CWID@med.cornell.edu
##SBATCH --mail-type=END

source /home/des2037/.bashrc

spack load -r /62q4vgx # this is cuda@9.2.88
#spack load -r cuda@9.2.88
spack load -r /omvpd5u # fftw@3.3.8

conda activate pHreplicaExchange 

pH=$1
iteration=$2
nsteps=$3
python -u run_replica.py $pH $iteration $nsteps
