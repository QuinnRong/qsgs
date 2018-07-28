gpu:
	sbatch tf_gpu.slurm

cpu:
	sbatch tf_cpu.slurm

fat:
	sbatch tf_fat.slurm

.PHONY:clean, delete_log, delete_model
clean:
	-rm *.err *.out
delete_log:
	-rm *.txt
delete_model:
	-rm -rf Model
