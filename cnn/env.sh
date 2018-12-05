source /usr/share/Modules/init/bash
module use /lustre/usr/modulefiles
module purge
module load miniconda3/4.3

source activate tf-py3-cpu
