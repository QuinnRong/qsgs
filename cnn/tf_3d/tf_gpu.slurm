#!/bin/bash

#SBATCH -J tf-2d-gpu
#SBATCH -p gpu
#SBATCH --mail-type=end
#SBATCH --mail-user=qyrong@sjtu.edu.cn
#SBATCH -o tf-%j.out
#SBATCH -e tf-%j.err
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --gres=gpu:2

source /usr/share/Modules/init/bash
module use /lustre/usr/modulefiles
module purge
module load miniconda3/4.3 cuda/9.0 cudnn/7.0

source activate tf-py3-gpu
python tf_3d.py
