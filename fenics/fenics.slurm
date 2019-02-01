#!/bin/bash

#SBATCH -J fen0
#SBATCH -p cpu
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

python poisson3d.py run_1_valid 0 0 100 100
python poisson3d.py run_1_valid 1 100 200 100
python poisson3d.py run_1_valid 2 200 300 100
python poisson3d.py run_1_valid 3 300 400 100
python poisson3d.py run_1_valid 4 400 500 100