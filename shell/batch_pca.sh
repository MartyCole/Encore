DIR=/nas/longleaf/home/mrcole/Encore/

sbatch --output="../logs/css_pca_cog_N050_reg.out" pca_reg.sh ${DIR}/results/abcd_subject_lists/cog_N050.mat ${DIR}/results/pca/reg_cog_N050 /work/users/m/r/mrcole/ScratchingTheSurface/new_reg/
sbatch --output="../logs/css_pca_cog_N100_reg.out" pca_reg.sh ${DIR}/results/abcd_subject_lists/cog_N100.mat ${DIR}/results/pca/reg_cog_N100 /work/users/m/r/mrcole/ScratchingTheSurface/new_reg/
sbatch --output="../logs/css_pca_cog_N200_reg.out" pca_reg.sh ${DIR}/results/abcd_subject_lists/cog_N200.mat ${DIR}/results/pca/reg_cog_N200 /work/users/m/r/mrcole/ScratchingTheSurface/new_reg/
sbatch --output="../logs/css_pca_cog_N400_reg.out" pca_reg.sh ${DIR}/results/abcd_subject_lists/cog_N400.mat ${DIR}/results/pca/reg_cog_N400 /work/users/m/r/mrcole/ScratchingTheSurface/new_reg/

sbatch --output="../logs/css_pca_cog_N050.out" pca.sh ${DIR}/results/abcd_subject_lists/cog_N050.mat ${DIR}/results/pca/cog_N050 
sbatch --output="../logs/css_pca_cog_N100.out" pca.sh ${DIR}/results/abcd_subject_lists/cog_N100.mat ${DIR}/results/pca/cog_N100 
sbatch --output="../logs/css_pca_cog_N200.out" pca.sh ${DIR}/results/abcd_subject_lists/cog_N200.mat ${DIR}/results/pca/cog_N200 
sbatch --output="../logs/css_pca_cog_N400.out" pca.sh ${DIR}/results/abcd_subject_lists/cog_N400.mat ${DIR}/results/pca/cog_N400 

sbatch --output="../logs/css_pca_cog_singlesex_N050_reg.out" pca_reg.sh ${DIR}/results/abcd_subject_lists/cog_singlesex_N050.mat ${DIR}/results/pca/reg_cog_singlesex_N050 /work/users/m/r/mrcole/ScratchingTheSurface/new_reg/
sbatch --output="../logs/css_pca_cog_singlesex_N100_reg.out" pca_reg.sh ${DIR}/results/abcd_subject_lists/cog_singlesex_N100.mat ${DIR}/results/pca/reg_cog_singlesex_N100 /work/users/m/r/mrcole/ScratchingTheSurface/new_reg/
sbatch --output="../logs/css_pca_cog_singlesex_N200_reg.out" pca_reg.sh ${DIR}/results/abcd_subject_lists/cog_singlesex_N200.mat ${DIR}/results/pca/reg_cog_singlesex_N200 /work/users/m/r/mrcole/ScratchingTheSurface/new_reg/
sbatch --output="../logs/css_pca_cog_singlesex_N400_reg.out" pca_reg.sh ${DIR}/results/abcd_subject_lists/cog_singlesex_N400.mat ${DIR}/results/pca/reg_cog_singlesex_N400 /work/users/m/r/mrcole/ScratchingTheSurface/new_reg/

sbatch --output="../logs/css_pca_cog_singlesex_N050.out" pca.sh ${DIR}/results/abcd_subject_lists/cog_singlesex_N050.mat ${DIR}/results/pca/cog_singlesex_N050 
sbatch --output="../logs/css_pca_cog_singlesex_N100.out" pca.sh ${DIR}/results/abcd_subject_lists/cog_singlesex_N100.mat ${DIR}/results/pca/cog_singlesex_N100 
sbatch --output="../logs/css_pca_cog_singlesex_N200.out" pca.sh ${DIR}/results/abcd_subject_lists/cog_singlesex_N200.mat ${DIR}/results/pca/cog_singlesex_N200 
sbatch --output="../logs/css_pca_cog_singlesex_N400.out" pca.sh ${DIR}/results/abcd_subject_lists/cog_singlesex_N400.mat ${DIR}/results/pca/cog_singlesex_N400 

sbatch --output="../logs/css_pca_sex_N050_reg.out" pca_reg.sh ${DIR}/results/abcd_subject_lists/sex_N050.mat ${DIR}/results/pca/reg_sex_N050 /work/users/m/r/mrcole/ScratchingTheSurface/new_reg/
sbatch --output="../logs/css_pca_sex_N100_reg.out" pca_reg.sh ${DIR}/results/abcd_subject_lists/sex_N100.mat ${DIR}/results/pca/reg_sex_N100 /work/users/m/r/mrcole/ScratchingTheSurface/new_reg/
sbatch --output="../logs/css_pca_sex_N200_reg.out" pca_reg.sh ${DIR}/results/abcd_subject_lists/sex_N200.mat ${DIR}/results/pca/reg_sex_N200 /work/users/m/r/mrcole/ScratchingTheSurface/new_reg/
sbatch --output="../logs/css_pca_sex_N400_reg.out" pca_reg.sh ${DIR}/results/abcd_subject_lists/sex_N400.mat ${DIR}/results/pca/reg_sex_N400 /work/users/m/r/mrcole/ScratchingTheSurface/new_reg/

sbatch --output="../logs/css_pca_sex_N050.out" pca.sh ${DIR}/results/abcd_subject_lists/sex_N050.mat ${DIR}/results/pca/sex_N050 
sbatch --output="../logs/css_pca_sex_N100.out" pca.sh ${DIR}/results/abcd_subject_lists/sex_N100.mat ${DIR}/results/pca/sex_N100 
sbatch --output="../logs/css_pca_sex_N200.out" pca.sh ${DIR}/results/abcd_subject_lists/sex_N200.mat ${DIR}/results/pca/sex_N200 
sbatch --output="../logs/css_pca_sex_N400.out" pca.sh ${DIR}/results/abcd_subject_lists/sex_N400.mat ${DIR}/results/pca/sex_N400 

