#!/bin/bash

#SBATCH --job-name=QSGS
#SBATCH --partition=cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=qyrong@sjtu.edu.cn
#SBATCH -n 1
#SBATCH --output=QSGS-%j.out
#SBATCH --error=QSGS-%j.err

ulimit -s unlimited         
ulimit -l unlimited 

source /usr/share/Modules/init/bash
module purge
module load icc impi

./qsgs_3d run_1_valid 500 100
./qsgs_3d run_2_train 2000 100
./qsgs_3d run_3_test 500 100
