% SCRIPT NAME:
%   single_subject_simulation
%
% DESCRIPTION:
%   Perform a test registration in both directions, Q1 -> Q2 and Q2 >- Q1. 
%
% MATLAB VERSION:
%   R2022b
%

restoredefaultpath
addpath(genpath('../../'))

clear all
close all

%% Setup

rng(23711228);                         % set seed for reproduciblity

ICO_RESOLUTION = 4;                    % resolution settings
l = 30;                                % order of SH for tangent basis

ico_mesh = icosphere(ICO_RESOLUTION);  % generate a triagulated mesh for the grid
grid = SphericalGrid(ico_mesh, l);     % generate a grid object from the mesh

%% Simulate connectivity on the icosphere

% load real data example (endpoints on a sphere)
tmp = load("../data/100206_native_coords.mat");

F1_hemi_in = tmp.surf_in;
F1_hemi_out = tmp.surf_out;

F1_start_pts = tmp.vtx_in;
F1_end_pts = tmp.vtx_out;

% the concon class object representing connectomes
F1 = Concon(grid,grid,F1_start_pts,F1_end_pts,F1_hemi_in,F1_hemi_out);

% load real data example (endpoints on a sphere)
tmp = load("../data/106824_native_coords.mat");

F2_hemi_in = tmp.surf_in;
F2_hemi_out = tmp.surf_out;

F2_start_pts = tmp.vtx_in;
F2_end_pts = tmp.vtx_out;

% the concon class object representing connectomes
F2 = Concon(grid,grid,F2_start_pts,F2_end_pts,F2_hemi_in,F2_hemi_out);

%% Create the smoothing Kernel

% calculate heat kernel for smoothing connectomes
kernel_builder = SphericalHeatKernel(grid,grid,l);

[sigma, ISE] = kernel_builder.cross_validate_sigma(F1,linspace(0.0005,0.005,10));

figure(50)
plot(linspace(0.0005,0.005,10),ISE)
ylabel('ISE')
xlabel('Sigma')
title(sprintf("Optimal Sigma: %0.4f",sigma))

[K, dK] = kernel_builder.compute_sigma(sigma);

%% Registration test 

% the registration algorithm class object (note that smaller l encourages
% smoother diffeomorphisms due to the simpler number of directions)
encore = Encore(grid,grid,15,0.05,100,1e-4);

% register in both directions
[F1_lh_warp,F1_rh_warp,~] = encore.register(F2,F1,K,dK,'init_rotation',true,'verbose',11);
[F2_lh_warp,F2_rh_warp,~] = encore.register(F1,F2,K,dK,'init_rotation',true,'verbose',12);

%% Results

orig_F1 = gather(F1.evaluate(K));
orig_F2 = gather(F2.evaluate(K));

F1_to_F2 = F1.copy();
F1_to_F2.warp_connectome(F1_lh_warp, F1_rh_warp);

F2_to_F1 = F2.copy();
F2_to_F1.warp_connectome(F2_lh_warp, F2_rh_warp);

reg_F1toF2 = gather(F1_to_F2.evaluate(K));
reg_F2toF1 = gather(F2_to_F1.evaluate(K));

A = [grid.A; grid.A] * [grid.A; grid.A]';

fprintf('Q2 -> Q1\nOriginal Distance: %0.3f\nF1->F2 Distance: %0.3f\nF2->F1 Distance: %0.2f\n\n', ...
     acos(sum(A.*sqrt(orig_F1).*sqrt(orig_F2),'all')), ...
     acos(sum(A.*sqrt(reg_F1toF2).*sqrt(orig_F2),'all')), ...
     acos(sum(A.*sqrt(reg_F2toF1).*sqrt(orig_F1),'all')))

fig = figure(13);
t = tiledlayout(1,3);
t.Padding = 'compact';
t.TileSpacing = 'compact';

nexttile
imagesc(log(A.*(orig_F1 - orig_F2).^2+1));
clim([0,5e-7])
xlim([4302,4438]);
ylim([4302,4438]);
xticklabels([]);
yticklabels([]);
title(sprintf('Original (F1 - F2)^2 Dist = %0.3f', acos(sum(A.*sqrt(orig_F1).*sqrt(orig_F2),'all'))));
nexttile
imagesc(log(A.*(reg_F1toF2 - orig_F2).^2+1))
clim([0,5e-7])
xlim([4302,4438]);
ylim([4302,4438]);
xticklabels([]);
yticklabels([]);
title(sprintf('Registered (F1 - F2)^2 Dist = %0.3f', acos(sum(A.*sqrt(reg_F1toF2).*sqrt(orig_F2),'all'))));
nexttile
imagesc(log(A.*(reg_F2toF1 - orig_F1).^2+1))
clim([0,5e-7])
xlim([4302,4438]);
ylim([4302,4438]);
xticklabels([]);
yticklabels([]);
colorbar()
title(sprintf('Registered (F2 - F1)^2 Dist = %0.3f', acos(sum(A.*sqrt(reg_F2toF1).*sqrt(orig_F1),'all'))));
