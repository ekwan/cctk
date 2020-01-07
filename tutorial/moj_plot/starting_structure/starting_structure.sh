#!/bin/bash 

#SBATCH -N 1        # number of nodes

#SBATCH -n 16   # number of cores 

#SBATCH -p defq   # partition to submit to 

#SBATCH --mem=34000  # memory per node in MB (see also --mem-per-cpu)

#SBATCH -t 1440 # expected runtime in minutes

#SBATCH -J g16_starting_structure # name of this job

# This is a template script for submitting Gaussian jobs
# to SLURM.  Tags starting with @ will be replaced with
# a script.  Eugene Kwan, May 2014

# prevent core dumps on job failure
ulimit -c 0

# set scratch
mkdir $LOCAL_SCRATCH/cwagen_starting_structure/
export GAUSS_SCRDIR=$LOCAL_SCRATCH/cwagen_starting_structure

# write out when and where the job started

echo "*************************************" > log.txt
echo "Scratch is: " $LOCAL_SCRATCH/cwagen_starting_structure >> log.txt
echo "Running on host:" >> log.txt
hostname >> log.txt
echo "Job starting_structure started at..." >> log.txt
date >> log.txt

# run job
g16 starting_structure.gjf starting_structure.out

# remove scratch
rm -rf $LOCAL_SCRATCH/cwagen_starting_structure

# analyze the result
./analyze.sh starting_structure.out >> log.txt

# add it to the master log
cat log.txt >> ../output/output.txt

# move completed job to output directory
mv starting_structure.out "../output/"


