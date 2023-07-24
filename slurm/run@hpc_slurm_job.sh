#!/bin/bash

##CHANGE THIS!

RUN_NAME="run-PROT4-bench_dataset@hpc"

##SBATCH --time=0:10:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

DATASET_NAME="dados"
INPUT_DATASET_FILE="../../datasets/bench_dataset.hd5.original"

## MAYBE CHANGE THIS!

EXE="../bin/laid-hdf5-mpi"

## DON'T CHANGE THIS!

# Be sure to request the correct partition to avoid the job to be held in the queue, furthermore
#	on CIRRUS-B (Minho)  choose for example HPC_4_Days
#	on CIRRUS-A (Lisbon) choose for example hpc
#SBATCH --partition=hpc

# Used to guarantee that the environment does not have any other loaded module
module purge

# Load software modules. Please check session software for the details
module load gcc11/libs/hdf5/1.14.0

# Disable warning for mismatched library versions
# Cirrus.8 has different hdf5 versions on short and hpc partitions
# even if we load the same module
# ##Headers are 1.14.0, library is 1.10.5
HDF5_DISABLE_VERSION_CHECK=1
export HDF5_DISABLE_VERSION_CHECK

# Run
echo "=== Running ==="
if [ -f "$EXE" ]; then
    chmod u+x $EXE
    mpiexec --display-map -np $SLURM_NTASKS $EXE -d $DATASET_NAME -f $INPUT_DATASET_FILE
else
    echo "$EXE Not found!"
fi

echo "Finished with job $SLURM_JOBID"
