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

N = 500000;       % number of fibers to generate
kappa1 = 1;       % concentration term (higher means less spread)
kappa2 = 2;       % concentration term (higher means less spread)
mu1 = [pi/2,0];   % mean point to create start points around
mu2 = [pi/2,pi];  % mean point to create end points around

% generate a random warp
lh_warp = random_diffeomorphism(grid,8,0.75,8000);
rh_warp = random_diffeomorphism(grid,8,0.75,8000);

% generate a random rotation
lh_rotation = random_rotation(3);
rh_rotation = random_rotation(3);

lh_warp_and_rot = lh_warp.copy();
rh_warp_and_rot = lh_warp.copy();
lh_warp_and_rot.V = lh_warp.V * lh_rotation;
rh_warp_and_rot.V = rh_warp.V * rh_rotation;

% generate the connectomes
[orig_sp,orig_ep] = random_connectome(N,kappa1,mu1,kappa2,mu2);

hemi1 = randsample([0 1], size(orig_ep,2), true);
hemi2 = randsample([0 1], size(orig_ep,2), true);
hemi2(1:floor(0.8*size(orig_ep,2))) = hemi1(1:floor(0.8*size(orig_ep,2)));

F1 = Concon(grid,grid,orig_sp.',orig_ep.',hemi1,hemi2);
F2 = copy(F1);
F2.warp_connectome(lh_warp,rh_warp);

F3 = copy(F1);
F3.warp_connectome(lh_warp_and_rot,rh_warp_and_rot);

%% Calculate the heat kernel

% calculate heat kernel for smoothing connectomes
kernel_builder = SphericalHeatKernel(grid,grid,l);
[K, dK] = kernel_builder.compute_sigma(0.005);

%% Registration test 

% the registration algorithm class object (note that smaller l encourages
% smoother diffeomorphisms due to the simpler number of directions)
encore = Encore(grid,grid,30,0.05,100,1e-8);

% register in both directions
[F1_lh_warp,F1_rh_warp,~] = encore.register(F2,F1,K,dK,'init_rotation',false,'verbose',1);
[F1_lh_warp_and_rot,F1_rh_warp_and_rot,~] = encore.register(F3,F1,K,dK,'init_rotation',true,'verbose',2);

%% Results

orig_F1 = gather(F1.evaluate(K));
orig_F2 = gather(F2.evaluate(K));

F1_warped = F1.copy();
F1_warped.warp_connectome(F1_lh_warp, F1_rh_warp);

reg_F1toF2 = gather(F1_warped.evaluate(K));

A = [grid.A; grid.A] * [grid.A; grid.A]';
fprintf('Q2 -> Q1\nOriginal Distance: %0.2f\nRegistered Distance: %0.2f\n\n', ...
     acos(sum(A.*sqrt(orig_F1).*sqrt(orig_F2),'all')), ...
     acos(sum(A.*sqrt(reg_F1toF2).*sqrt(orig_F2),'all')))
 
fig = figure(3);
t = tiledlayout(1,2);
t.Padding = 'compact';
t.TileSpacing = 'compact';

nexttile
lh_warp.plot('True LH-Warp')
clim([0 2])
view([0 90 0]);
nexttile
rh_warp.plot('True RH-Warp')
clim([0 2])
view([0 90 0]);

fig.Units = "normalized";
fig.OuterPosition = [0.154166666666667 0.08 0.411805555555556 0.634444444444444];

%%

fig = figure(4);
t = tiledlayout(1,2);
t.Padding = 'compact';
t.TileSpacing = 'compact';

nexttile
imagesc(((orig_F1 + 1) - (orig_F2 + 1)).^2);
clim([0,5e-5])
xlim([4600,5124]);
ylim([4600,5124]);
xticklabels([]);
yticklabels([]);
title(sprintf('Original (F1 - F2)^2 Dist = %0.2f', acos(sum(A.*sqrt(orig_F1).*sqrt(orig_F2),'all'))));
nexttile
imagesc(((reg_F1toF2 + 1) - (orig_F2 + 1)).^2)
clim([0,5e-5])
xlim([4600,5124]);
ylim([4600,5124]);
xticklabels([]);
yticklabels([]);
colorbar()
title(sprintf('Registered (F1 - F2)^2 Dist = %0.2f', acos(sum(A.*sqrt(reg_F1toF2).*sqrt(orig_F2),'all'))));

%%
orig_F3 = gather(F3.evaluate(K));

F3_warped = F1.copy();
F3_warped.warp_connectome(F1_lh_warp_and_rot, F1_rh_warp_and_rot);

reg_F1toF3 = gather(F3_warped.evaluate(K));

fig = figure(5);
t = tiledlayout(1,2);
t.Padding = 'compact';
t.TileSpacing = 'compact';

nexttile
imagesc(((orig_F1 + 1) - (orig_F3 + 1)).^2);
clim([0,5e-5])
xlim([4600,5124]);
ylim([4600,5124]);
xticklabels([]);
yticklabels([]);
title(sprintf('Original (F1 - F3)^2 Dist = %0.2f', acos(sum(A.*sqrt(orig_F1).*sqrt(orig_F3),'all'))));
nexttile
imagesc(((reg_F1toF3 + 1) - (orig_F3 + 1)).^2)
clim([0,5e-5])
xlim([4600,5124]);
ylim([4600,5124]);
xticklabels([]);
yticklabels([]);
colorbar()
title(sprintf('Registered (F1 - F3)^2 Dist = %0.2f', acos(sum(A.*sqrt(reg_F1toF3).*sqrt(orig_F3),'all'))));