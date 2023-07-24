%
%  File    :   register_subjects.m 
%  Author  :   Martin R. Cole
%  Date    :   19/07/2023
%  

function register_subjects(pFrom, pTo, pOutput)

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
input_file = "mesh_intersections_ico4.mat";
output_folder = pOutput;

%% Compute registration

lh_mesh.V = normr(ico_grid.lh_V);
lh_mesh.T = ico_grid.lh_T;
rh_mesh.V = normr(ico_grid.rh_V);
rh_mesh.T = ico_grid.rh_T;

encore = Encore(lh_mesh,rh_mesh,30,0.01,1000,0);

for i = pFrom:pTo   
    fprintf("Registration for %0.4i\n",i);

    % load the oversmoothed template
    load('/work/users/m/r/mrcole/ScratchingTheSurface/new_reg/template_N200_FWHM40_raw.mat','template');

    % load the subject endpoint data
    tmp = load(sprintf('%s/%s/%s/%s', input_folder, selected_subs(i), input_subfolder, input_file));
    
    % get the cartesian coordinates of the endpoints (currently saved as
    % barycentric coordinates, but I will probably change that so we can
    % skip this step in a more finalised version of the pipeline
    [coords_in, coords_out, surf_in, surf_out] = get_endpoint_coords(ico_grid, tmp);
    coords_in = normr(coords_in);
    coords_out = normr(coords_out);

    % initialise the concon object (overly smooth kernel)
    concon = Concon(encore.lh_grid,encore.rh_grid,[coords_in,surf_in'],[coords_out,surf_out'],40,0.01,1e-10);
 
    fprintf("=============================\n");
    fprintf("Starting initial registration\n");
    fprintf("=============================\n");
    tic
    % perform initial registration to get a smooth global warp
    [lh_warp,rh_warp,cost1] = encore.register(template,concon,true);
    toc  

    % load the final template
    load('/work/users/m/r/mrcole/ScratchingTheSurface/new_reg/template_N200_FWHM20_raw.mat','template');

    % set the kernel to the value we will use in the final analysis
    concon = concon.set_kernel(20, 0.01);

    fprintf("===============================\n");
    fprintf("Starting secondary registration\n");
    fprintf("===============================\n");
    tic
    % perform second registration (starting from previous results) to get a smooth local warp
    [lh_warp,rh_warp,cost2] = encore.register(template,concon,true,lh_warp,rh_warp);
    toc    

    fprintf("===============================\n");
    fprintf("Finished registration for %0.4i\n",i);
    fprintf("Cost1: %0.4f\nCost2: %0.4f\n",cost1,cost2);
    fprintf("===============================\n\n");

    save(sprintf('%s/css_registered_warp%0.4i.mat', output_folder, i), 'lh_warp', 'rh_warp', '-v7.3'); 
end

end