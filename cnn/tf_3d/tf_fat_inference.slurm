#!/bin/bash

#SBATCH -J tf-2d-fat
#SBATCH -p fat
#SBATCH --mail-type=end
#SBATCH --mail-user=qyrong@sjtu.edu.cn
#SBATCH -o tf-%j.out
#SBATCH -e tf-%j.err
#SBATCH -n 1
#SBATCH --exclusive

source /usr/share/Modules/init/bash
module use /lustre/usr/modulefiles
module purge
module load miniconda3/4.3

source activate tf-py3-cpu
python tf_3d_inference.py
