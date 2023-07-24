#!/bin/bash
#SBATCH --time=0-08:00:00
#SBATCH --mem=180GB

module load matlab/2022b
cd ..
matlab -nodisplay -singleCompThread -r "get_template; exit" 
