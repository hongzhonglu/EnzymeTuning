#!/bin/bash

#SBATCH --job-name=gan_test
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

cd /dssg/home/acct-clslhz/clslhz/WXT_EC/tensorflow


python Step2_gan_training2.py

#python test_gan.py