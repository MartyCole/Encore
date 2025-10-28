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

N = 500000;       % number of fibers to generate
kappa1 = 1;       % concentration term (higher means less spread)
kappa2 = 2;       % concentration term (higher means less spread)
mu1 = [pi/2,0];   % mean point to create start points around
mu2 = [pi/2,pi];  % mean point to create end points around

% generate a random warp
lh_warp = random_diffeomorphism(grid,8,0.75,8000);
rh_warp = random_diffeomorphism(grid,8,0.75,8000);

% generate the connectomes
[orig_sp,orig_ep] = random_connectome(N,kappa1,mu1,kappa2,mu2);

hemi1 = randsample([0 1], size(orig_ep,2), true);
hemi2 = randsample([0 1], size(orig_ep,2), true);
hemi2(1:(0.8*size(orig_ep,2))) = hemi1(1:(0.8*size(orig_ep,2)));

F1 = Concon(grid,grid,orig_sp.',orig_ep.',hemi1,hemi2);
F2 = copy(F1);
F2.warp_connectome(lh_warp,rh_warp);

%% Calculate the heat kernel

% calculate smoothed connectomes
kernel = spherical_heat_kernel(grid,grid,0.005,l);

%% Registration test 

encore = Encore(grid,grid,l,0.05,1000,1e-10);

[Freg,est_lh_warp,est_rh_warp,cost] = encore.register(F2,F1,kernel,'verbose',1);

%% Results

orig_F1 = F1.evaluate(kernel, true);
orig_F2 = F2.evaluate(kernel, true);
reg_F1toF2 = Freg.evaluate(kernel, true);

fprintf('Q2 -> Q1\nOriginal Distance: %0.2f\nRegistered Distance: %0.2f\n\n', ...
     acos(sum(sqrt(orig_F1).*sqrt(orig_F2),'all')), ...
     acos(sum(sqrt(reg_F1toF2).*sqrt(orig_F2),'all')))
 
fig = figure(2);
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

fig = figure(3);
t = tiledlayout(1,2);
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
colorbar()
title(sprintf('Registered (F1 - F2)^2 Dist = %0.2f', acos(sum(sqrt(reg_F1toF2).*sqrt(orig_F2),'all'))));