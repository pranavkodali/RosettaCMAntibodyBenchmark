#!/bin/bash

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=pranav.kodali@vanderbilt.edu
#SBATCH --job-name=PATHL
#SBATCH --output=PATHL_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=10-0
#SBATCH --array=0-99%THREADS
#exclude=bf21,bf22,bf23,bf24

//export LD_LIBRARY_PATH="/dors/meilerlab/apps/Linux2/x86_64/lib64:/dors/meilerlab/apps/Linux2/x86_64/lib"

echo "start"
echo ${SLURM_ARRAY_TASK_ID}
echo "end"
bash ./sc_parser.bash "${SLURM_ARRAY_TASK_ID}_score.sc"

module load GCC/5.4.0-2.26

/dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin/relax.default.linuxgccrelease @OPTIONSP -l "${SLURM_ARRAY_TASK_ID}_score.txt" -out:prefix relax -nstruct 1
