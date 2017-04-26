This is the final Project For MATH 440


For MIO:
module purge
module load openmpi/gcc/default
make
sbatch -p mganesh test.slurm
squeue -u<yourusername>
