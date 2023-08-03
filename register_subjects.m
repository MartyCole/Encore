%
%  File    :   register_subjects.m 
%  Author  :   Martin R. Cole
%  Date    :   19/07/2023
%  

function register_subjects(pFrom, pTo, pOutput)

addpath('interpolation/')
addpath('kde/')
addpath('simulations/')
addpath('tangent_basis')
addpath('utils/')

%% Global settings

% grid to do evaluations over
grid_mesh = load('/nas/longleaf/home/mrcole/ScratchingTheSurface/Code/SBCI_AVG/template_sphere_grid_ico4.mat');

% selected_subs are the IDs of the subjects we're interested in
load('/nas/longleaf/home/mrcole/ScratchingTheSurface/Code/exampleData/selected_subs_1000_totalcog.mat')

% location of the connectivity data we would like to load
input_folder = "/overflow/zzhanglab/SBCI_Finished_ABCD_Data/";
input_subfolder = "/psc_sbci_final_files/sbci_connectome/";
input_file = "smoothed_sc_avg_0.005_ico4.mat";
output_folder = pOutput;

load(sprintf('/work/users/m/r/mrcole/ScratchingTheSurface/new_reg/template_N200_control_log.mat'));

%% Compute registration

lh_mesh.V = normr(grid_mesh.lh_V);
lh_mesh.T = grid_mesh.lh_T;
rh_mesh.V = normr(grid_mesh.rh_V);
rh_mesh.T = grid_mesh.rh_T;

encore = Encore(lh_mesh,rh_mesh,30,0.01,1000,0);

for i = pFrom:pTo
    tmp = load(sprintf('%s/%s/%s/%s', input_folder, selected_subs(i), input_subfolder, input_file));

    Q = (tmp.sc + tmp.sc');
    Q = Q ./ sum(tmp.sc(:));    
    Q = log(1e6*Q + 1); 

    tic
    [SC,lh_warp,rh_warp,~] = encore.register(template,Q,true);
    toc
    
    save(sprintf('%s/registered_sub%0.4i.mat', output_folder, i), 'SC', 'lh_warp', 'rh_warp', '-v7.3'); 
end

end
