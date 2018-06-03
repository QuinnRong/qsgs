#!/bin/bash

#SBATCH --job-name=QSGS
#SBATCH --partition=fat
#SBATCH --mail-type=end
#SBATCH --mail-user=qyrong@sjtu.edu.cn
#SBATCH -n 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

ulimit -s unlimited         
ulimit -l unlimited 

source /usr/share/Modules/init/bash
module use /lustre/usr/modulefiles/pi
module purge
module load icc impi

./a.out
