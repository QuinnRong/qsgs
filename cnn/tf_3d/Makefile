gpu:
	sbatch tf_gpu.slurm

cpu:
	sbatch tf_cpu.slurm

fat:
	sbatch tf_fat.slurm

test_cpu:
	sbatch tf_cpu_inference.slurm

test_gpu:
	sbatch tf_gpu_inference.slurm

test_fat:
	sbatch tf_fat_inference.slurm

.PHONY:clean, clean_all
clean:
	-rm *.err *.out
clean_all:
	-rm *.err *.out *.txt
	-rm -rf Model
