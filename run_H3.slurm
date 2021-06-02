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

export LD_LIBRARY_PATH=/dors/meilerlab/apps/Linux2/x86_64/pyrosetta/PyRosetta.ScientificLinux-r56324.64Bit/rosetta:/dors/meilerlab/apps/Linux2/x86_64/lib64:/dors/meilerlab/apps/Linux2/x86_64/lib:/sb/apps/Linux/x86_64/lib:/sb/apps/Linux/lib:/dors/meilerlab/apps/Linux2/x86_64/python/external_modules_src/pymol-centos-1.7.0.4/lib64:/dors/meilerlab/apps/Linux2/x86_64/python/external_modules_src/pymol-centos-1.7.0.4/lib:/sb/apps/Linux/x86_64/lib:/sb/apps/Linux/lib:/dors/meilerlab/apps/Linux2/x86_64/gcc/5.2.0/lib64/:/home/kodalip/Pranav/RosettaTutorials/Rosetta_AB/main/source/build/external/release/linux/3.10/64/x86/gcc/5.2/default/:/home/kodalip/Pranav/RosettaTutorials/Rosetta_AB/main/source/build/external/release/linux/3.10/64/x86/gcc/5.2/default/:/dors/meilerlab/apps/Linux2/x86_64/lib64/:$LD_LIBRARY_PATH

module load GCC/5.4.0-2.26

/dors/meilerlab/home/kodalip/Pranav/RosettaTutorials/Rosetta_AB/main/source/bin/antibody_H3.linuxgccrelease @/dors/meilerlabOPTIONSP/abH3.flags -in:file:s /dors/meilerlabOPTIONSP/grafting/model-0.relaxed.pdb -nstruct 10 -antibody:auto_generate_kink_constraint -antibody:all_atom_mode_kink_constraint -out:file:scorefile /dors/meilerlabOPTIONSP/H3_modeling_scores.fasc -out:path:pdb /dors/meilerlabOPTIONSP/H3_modeling/ -out:prefix ${SLURM_ARRAY_TASK_ID}_
