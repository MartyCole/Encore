% SCRIPT NAME:
%   example_analyis
%
% DESCRIPTION:
%   Sample analysis using 10 participants of the HCP Young Adults dataset
%
% MATLAB VERSION:
%   R2022b
%

restoredefaultpath
addpath(genpath('../../'))

clear variables; close all

%% Setup

% set seed for reproduciblity
rng(23711228);                         

% resolution setting
ICO_RESOLUTION = 4;                   

% order of spherical harmonics for tangent basis
l = 30;

% a place to save the results
if ~exist("../example_results/","dir")
    mkdir("../example_results");
end

%% Setup the Analysis

fprintf("-------------------------------------------\n");
fprintf("Generating grid for registration\n")

% generate a triagulated mesh for the grid
ico_mesh = icosphere(ICO_RESOLUTION);  

% the grid to register with (both hemis)
grid = SphericalGrid(ico_mesh, l);    

% load subject data
sublist = readlines("../data/sublist", "EmptyLineRule", "skip");
N = length(sublist);
Fs = cell(N,1);

fprintf("Loading subject data:\n")
for i = 1:N
    fprintf("Sub: %d/%d\n",i,N)

    % load real data example (endpoints on a sphere)
    tmp = load(sprintf("../data/%s_reg_sphere_intersections.mat",sublist(i)));

    F1_hemi_in = double(tmp.surf_in);
    F1_hemi_out = double(tmp.surf_out);

    F1_start_pts = double(tmp.vtx_in);
    F1_end_pts = double(tmp.vtx_out);

    % the concon class object representing connectomes
    Fs{i} = SConcon(grid,grid,F1_start_pts,F1_end_pts,F1_hemi_in,F1_hemi_out);
end

%% Create the Grid and smoothing Kernel

fprintf("-------------------------------------------\n");
fprintf("Generating kernel for registration\n")

% calculate heat kernel for smoothing connectomes
kernel_builder = SphericalHeatKernel(grid,grid,l);

% estimated optimal diffusion time (changes based on subject, 0.005 is good)
[sigma, LOO] = kernel_builder.cross_validate_sigma(Fs{2},linspace(0.001,0.008,22));

fprintf("Evalutating the kernel and its derivative\n")

% evalutate the kernel and its derivative with the given sigma
[K, dK] = kernel_builder.compute_sigma(0.005);

% save for future evaluation of concons
save("../example_results/ico4_heatkernel.mat","K","dK")

%% Perform the registration (including estimating a template to register to)

fprintf("-------------------------------------------\n");
fprintf("Generating a template from the population\n");

% the registration algorithm class object (note that smaller l encourages
% smoother diffeomorphisms due to the simpler number of directions)
encore = Encore(grid,grid,15,0.05,100,1e-4);

% create a template
template = encore.get_template(Fs,K,100);

% register the population to the template
for i = 1:N
    fprintf("Subject: %i/%i\n", i, N)
    [lh_warp,rh_warp,~] = encore.register(template,Fs{i},K,dK,"verbose",1);
    save(sprintf("../example_results/%s_warp.mat",sublist(i)),"lh_warp","rh_warp")
end

%% Show resulting warps

figure(2)
tiledlayout(1,2)
nexttile
lh_warp.plot_field("LH Warp")
pbaspect([1 1 1]); view(0, 0);
nexttile
rh_warp.plot_field("RH Warp")
pbaspect([1 1 1]); view(0, 0);