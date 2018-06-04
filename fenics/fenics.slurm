#!/bin/bash

#SBATCH -J fenics
#SBATCH -p fat
#SBATCH --mail-type=end
#SBATCH --mail-user=qyrong@sjtu.edu.cn
#SBATCH -o fenics-%j.out
#SBATCH -e fenics-%j.err
#SBATCH -n 1
#SBATCH --exclusive

source /usr/share/Modules/init/bash
module purge
module use /lustre/usr/modulefiles
module load miniconda2/4.3

source activate fenicsproject
python Poisson3D_Mostafa_quinn.py
