#!/bin/bash -x
#SBATCH --job-name="test"
#comment	="Hello..."
#SBATCH --nodes=02
#SBATCH --ntasks-per-node=12
#SBATCH --ntasks=24
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=00:59:59
#SBATCH -o project.out
#SBATCH -e project.err

# Go to the directoy from which our job was launched
cd $SLURM_SUBMIT_DIR
# Create a short JOBID base on the one provided by the scheduler
JOBID=‘echo $SLURM_JOBID’
# Save a copy of our environment and script
cat $0 > script.$JOBID
printenv  > env.$JOBID
# Run the job.
# The echo will go into the standard output for this job
# The standard output file will end up in the directory
# from which the job was launched.
srun par2_exe
