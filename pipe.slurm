#!/bin/bash
# HEADER for Parallel job using 80 processors:
#SBATCH --nodes=2           # number of nodes
#SBATCH --ntasks-per-node=10    # number of processors per node
#SBATCH --cpus-per-task=4       # number of cpus per task
#SBATCH -t 1:00:00      # run for 1 hr max
#SBATCH --mail-type=begin   # send email when process begins...
#SBATCH --mail-type=end     # ...and when it ends...
#SBATCH --mail-type=fail    # ...or when it fails.
#SBATCH --mail-user=your_email@princeton.edu # send notifications to this email
#SBATCH -e job.err              # Name of output file for error messages
#SBATCH -o job.out              # Name of output file for standard output

# BODY - commands to be run
# Load required modules
module load intel               
module load intel-mpi

# Set number of openmp threads and number of MPI tasks
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASKS=$(echo "$SLURM_NNODES*$(echo $SLURM_TASKS_PER_NODE | cut -d '(' -f 1)" | bc -l)

echo $OMP_NUM_THREADS
echo $NTASKS

# Print some information
echo Directory is `pwd`
echo Time is `date`
echo
echo This job runs on the following processors:
echo $SLURM_JOB_NODELIST
echo This job has allocated $SLURM_NNODES nodes with $SLURM_TASKS_PER_NODE cores per node.
echo

# run arts
srun -n $NTASKS ~/NGA/bin/arts input_pipe