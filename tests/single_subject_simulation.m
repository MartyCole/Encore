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

N = 1000000;      % number of fibers to generate
kappa1 = 1;       % concentration term (higher means less spread)
kappa2 = 2;       % concentration term (higher means less spread)
mu1 = [pi/2,0];   % mean point to create start points around
mu2 = [pi/2,pi];  % mean point to create end points around

% generate a random warp
warp = random_diffeomorphism(grid,8,0.75,8000);

% generate the connectomes
[orig_sp,orig_ep] = random_connectome(N,kappa1,mu1,kappa2,mu2);
[warp_sp,warp_ep] = warp_connectome(grid,warp,orig_sp,orig_ep);

hemi1 = randsample([0 1], size(orig_ep,2), true);
hemi2 = randsample([0 1], size(orig_ep,2), true);
hemi2(1:(0.8*size(orig_ep,2))) = hemi1(1:(0.8*size(orig_ep,2)));

% calculate smoothed connectomes
kernel = spherical_kernel(grid,10,0.1);

F1 = spherical_kde(grid, kernel, orig_sp.', orig_ep.', hemi1, hemi2);
F2 = spherical_kde(grid, kernel, warp_sp.', warp_ep.', hemi1, hemi2);

%% Registration test 

encore = Encore(grid,grid,30,0.01,1000,0);

[SC_reg,lh_warp,rh_warp,~] = encore.register(F1,F2,'verbose',1);

fprintf('Q2 -> Q1\nOriginal Distance: %0.2f\nRegistered Distance: %0.2f\n\n', ...
     sqrt(sum(diag((F1 - F2)*(F1 - F2)'))), ...
     sqrt(sum(diag((F1 - SC_reg)*(F1 - SC_reg)'))))
 
figure(1)
lh_warp.plot('Estimated LH-Warp')
clim([0 2])
view([0 90 0]);
figure(2)
warp.plot('True LH-Warp')
clim([0 2])
view([0 90 0]);

