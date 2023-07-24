#!/bin/bash
#SBATCH --time=0-20:00:00
#SBATCH --mem=10GB

module load matlab/2022b
cd ..
matlab -nodisplay -r "register_subjects(${1},${2},'/work/users/m/r/mrcole/ScratchingTheSurface/new_reg/'); exit" 
