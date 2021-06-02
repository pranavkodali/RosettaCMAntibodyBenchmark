#!/bin/bash

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=pranav.kodali@vanderbilt.edu
#SBATCH --job-name=PATHL
#SBATCH --output=PATHL_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1-15
#SBATCH --array=0-99%THREADS
#exclude=bf21,bf22,bf23,bf24

//export LD_LIBRARY_PATH="/dors/meilerlab/apps/Linux2/x86_64/lib64:/dors/meilerlab/apps/Linux2/x86_64/lib"

module load GCC/5.4.0-2.26

/dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin/rosetta_scripts.default.linuxgccrelease @/dors/meilerlabOPTIONSP -nstruct NUMBEROFTRIALS -out:prefix ${SLURM_ARRAY_TASK_ID}_
