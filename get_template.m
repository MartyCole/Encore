%
%  File    :   find_template.m 
%  Author  :   Martin R. Cole
%  Date    :   20/07/2023
%  

addpath('interpolation/')
addpath('simulations/')
addpath('tangent_basis')
addpath('utils/')
addpath('/nas/longleaf/home/mrcole/SBCI_Toolkit/io')

%% Global settings

% grid to do evaluations over
ico_grid = load('/nas/longleaf/home/mrcole/ABCD/SBCI_AVG/template_sphere_grid_ico4.mat');
% selected_subs are the IDs of the subjects we're interested in
load('/nas/longleaf/home/mrcole/ScratchingTheSurface/Code/exampleData/selected_subs_1000_totalcog.mat')

% location of the connectivity data we would like to load
input_folder = "/overflow/zzhanglab/SBCI_Finished_ABCD_Data/";
input_subfolder = "/psc_sbci_final_files/sbci_connectome/";
input_file = "smoothed_sc_avg_0.005_ico4.mat";

%% Find the template

% use the freesurfer fsaverage mesh as the grid
lh_mesh.V = normr(ico_grid.lh_V);
lh_mesh.T = ico_grid.lh_T;
rh_mesh.V = normr(ico_grid.rh_V);
rh_mesh.T = ico_grid.rh_T;

encore = Encore(lh_mesh,rh_mesh,30,0.01,1000,0);

N_subs = 200;

rng(54569783);
sub_idx = randi(1000,N_subs);

P = length(lh_mesh.V) + length(rh_mesh.V);
SC = zeros(P,P,N_subs);

for i = 1:N_subs    
    tmp = load(sprintf('%s/%s/%s/%s', input_folder, selected_subs(sub_idx(i)), input_subfolder, input_file));
  
    SC(:,:,i) = (tmp.sc + tmp.sc');
    SC(:,:,i) = SC(:,:,i) ./ sum(tmp.sc(:));    
    SC(:,:,i) = log(1e6*SC(:,:,i) + 1);

    disp(i);
end

template = encore.get_template(SC,200);

save("/work/users/m/r/mrcole/ScratchingTheSurface/new_reg/template_N200_log.mat","template",'-v7.3')

SC = zeros(P,P,N_subs);

for i = 1:N_subs    
    tmp = load(sprintf('%s/%s/%s/%s', input_folder, selected_subs(i+1), input_subfolder, input_file));
  
    SC(:,:,i) = (tmp.sc + tmp.sc');
    SC(:,:,i) = SC(:,:,i) ./ sum(tmp.sc(:));   
    SC(:,:,i) = log(1e6*SC(:,:,i) + 1);

    disp(i);
end

template = encore.get_template(SC,200);

save("/work/users/m/r/mrcole/ScratchingTheSurface/new_reg/template_N200_control_log.mat","template",'-v7.3')

sub_idx = [2:101,901:1000];
SC = zeros(P,P,N_subs);

for i = 1:N_subs    
    tmp = load(sprintf('%s/%s/%s/%s', input_folder, selected_subs(sub_idx(i)), input_subfolder, input_file));
  
    SC(:,:,i) = (tmp.sc + tmp.sc');
    SC(:,:,i) = SC(:,:,i) ./ sum(tmp.sc(:));    
    SC(:,:,i) = log(1e6*SC(:,:,i) + 1);

    disp(i);
end

template = encore.get_template(SC,200);

save("/work/users/m/r/mrcole/ScratchingTheSurface/new_reg/template_N200_both_log.mat","template",'-v7.3')
