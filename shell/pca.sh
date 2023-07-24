#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --mem=100GB

module load matlab/2022b
cd ../analysis/
matlab -nodisplay -singleCompThread -r "pca('${1}','${2}')" 
