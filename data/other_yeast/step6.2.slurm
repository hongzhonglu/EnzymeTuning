#!/bin/bash

#SBATCH --job-name=error_num
#SBATCH --partition=a100
#SBATCH -N 1
#SBATCH --ntasks-per-node=16

#SBATCH --gres=gpu:1          # use 1 GPU
#SBATCH --mail-type=end
#SBATCH --mail-user=ariel926@163.com
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate tf-env

cd /dssg/home/acct-clslhz/clslhz/WXT_EC


python Step4_gan_all_phen.py
