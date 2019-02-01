source activate fenicsproject

python poisson3d.py run_1_valid 0 0 2 20

python poisson3d.py run_2_train 0 0 2 20
python poisson3d.py run_2_train 1 2 4 20
python poisson3d.py run_2_train 2 4 6 20
python poisson3d.py run_2_train 3 6 8 20

python poisson3d.py run_3_test 0 0 2 20