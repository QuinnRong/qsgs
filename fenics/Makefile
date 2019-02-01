run:
	sbatch fenics.slurm

run_mac:
	bash mac.sh

combine:
	python combine.py run_1_valid 1
	python combine.py run_2_train 4
	python combine.py run_3_test 1

.PHONEY: clean
clean:
	-rm -rf *.err *.out output
