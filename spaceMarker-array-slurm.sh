#!/bin/bash
#SBATCH --job-name=SpaceMarker-Batch
#SBATCH --time=24:00:00
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --mail-type=end
#SBATCH --mail-user=szhan121@jhu.edu

#module load seurat/4.1.1
#module list

#### process job list data
filename=sample_list.dat
declare -a myArray
mapfile -t myArray < $filename
nrJobs=${#myArray[@]}

#IS_MARCC_JOB=false
IS_MARCC_JOB=true

R_SCRIPT="./spQSP_SplnMarkers_MARCC.R"

R CMD BATCH $R_SCRIPT ${myArray[$((SLURM_ARRAY_TASK_ID-1))]}

echo "Finished with job $SLURM_ARRAY_TASK_ID"
