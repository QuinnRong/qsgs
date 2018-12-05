source /usr/share/Modules/init/bash
module purge
module use /lustre/usr/modulefiles
module load miniconda2/4.3

source activate fenicsproject
