% SCRIPT NAME:
%   single_subject_simulation
%
% DESCRIPTION:
%   Perform a test registration in both directions, Q1 -> Q2 and Q2 >- Q1. 
%
% MATLAB VERSION:
%   R2022b
%

addpath('../')
addpath('../interpolation/')
addpath('../kde/')
addpath('../simulations/')
addpath('../tangent_basis')
addpath('../utils/')

%% Setup

rng(23711228);                         % set seed for reproduciblity

ICO_RESOLUTION = 4;                    % resolution settings
l = 30;                                % order of SH for tangent basis

ico_mesh = icosphere(ICO_RESOLUTION);  % generate a triagulated mesh for the grid
grid = SphericalGrid(ico_mesh, l);     % generate a grid object from the mesh

%% Simulate connectivity on the icosphere

% load real data example (endpoints on a sphere)
xing = readmatrix("./data/100206_xing_sphere_avg_coords.tsv",'FileType','text');

F1_hemi_in = 1-xing(:,2);
F1_hemi_out = 1-xing(:,7);

F1_start_pts = xing(:,3:5);
F1_end_pts = xing(:,8:10);

% the concon class object representing connectomes
F1 = Concon(grid,grid,F1_start_pts,F1_end_pts,F1_hemi_in,F1_hemi_out);

% load real data example (endpoints on a sphere)
xing = readmatrix("./data/106824_xing_sphere_avg_coords.tsv",'FileType','text');

F2_hemi_in = 1-xing(:,2);
F2_hemi_out = 1-xing(:,7);

F2_start_pts = xing(:,3:5);
F2_end_pts = xing(:,8:10);

% the concon class object representing connectomes
F2 = Concon(grid,grid,F2_start_pts,F2_end_pts,F2_hemi_in,F2_hemi_out);

% calculate heat kernel for smoothing connectomes
kernel = spherical_heat_kernel(grid,grid,0.005,l);

%% Registration test 

% the registration algorithm class object (note that smaller l encourages
% smoother diffeomorphisms due to the simpler number of directions)
encore = Encore(grid,grid,15,0.05,1000,1e-4);

% register in both directions
[F1_to_F2,F1_lh_warp,F1_rh_warp,~] = encore.register(F2,F1,kernel,'verbose',11);
[F2_to_F1,F2_lh_warp,F2_rh_warp,~] = encore.register(F1,F2,kernel,'verbose',12);

%% Results

orig_F1 = F1.evaluate(kernel, true);
orig_F2 = F2.evaluate(kernel, true);
reg_F1toF2 = F1_to_F2.evaluate(kernel, true);
reg_F2toF1 = F2_to_F1.evaluate(kernel, true);

fprintf('Q2 -> Q1\nOriginal Distance: %0.2f\nF1->F2 Distance: %0.2f\nF2->F1 Distance: %0.2f\n\n', ...
     acos(sum(sqrt(orig_F1).*sqrt(orig_F2),'all')), ...
     acos(sum(sqrt(reg_F1toF2).*sqrt(orig_F2),'all')), ...
     acos(sum(sqrt(reg_F2toF1).*sqrt(orig_F1),'all')))

fig = figure(13);
t = tiledlayout(1,3);
t.Padding = 'compact';
t.TileSpacing = 'compact';

nexttile
imagesc(((1e6*orig_F1 + 1) - (1e6*orig_F2 + 1)).^2);
clim([0,1])
xlim([4600,5124]);
ylim([4600,5124]);
xticklabels([]);
yticklabels([]);
title(sprintf('Original (F1 - F2)^2 Dist = %0.2f', acos(sum(sqrt(orig_F1).*sqrt(orig_F2),'all'))));
nexttile
imagesc(((1e6*reg_F1toF2 + 1) - (1e6*orig_F2 + 1)).^2)
clim([0,1])
xlim([4600,5124]);
ylim([4600,5124]);
xticklabels([]);
yticklabels([]);
title(sprintf('Registered (F1 - F2)^2 Dist = %0.2f', acos(sum(sqrt(reg_F1toF2).*sqrt(orig_F2),'all'))));
nexttile
imagesc(((1e6*reg_F2toF1 + 1) - (1e6*orig_F1 + 1)).^2)
clim([0,1])
xlim([4600,5124]);
ylim([4600,5124]);
xticklabels([]);
yticklabels([]);
colorbar()
title(sprintf('Registered (F2 - F1)^2 Dist = %0.2f', acos(sum(sqrt(reg_F2toF1).*sqrt(orig_F1),'all'))));