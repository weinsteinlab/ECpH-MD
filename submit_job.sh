#!/bin/sh
#SBATCH --job-name=RNase-pH
#SBATCH -p edison
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --mem=40G
##SBATCH --dependency=afterany:1142811
#SBATCH --mail-user=kots.katya@gmail.com
##SBATCH --mail-type=END

source /home/des2037/.bashrc

spack load -r cuda@9.2.88
spack load -r /omvpd5u # fftw@3.3.8

conda activate pHreplicaExchange 

#python -u exchange_min_replica.py
python -u Exchange-min-replica.py
