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
addpath('../graphics/')
addpath('../interpolation/')
addpath('../kde/')
addpath('../simulations/')
addpath('../tangent_basis')
addpath('../utils/')
addpath('../../export_fig/')

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
warp = random_diffeomorphism(grid,8,0.75,8000);

% generate the connectomes
[orig_sp,orig_ep] = random_connectome(N,kappa1,mu1,kappa2,mu2);
[warp_sp,warp_ep] = warp_connectome(grid,warp,orig_sp,orig_ep);

% calculate smoothed connectomes
kernel = spherical_kernel(grid,10,0.1);
Q1 = spherical_kde(grid, kernel, orig_sp.', orig_ep.');
Q2 = spherical_kde(grid, kernel, warp_sp.', warp_ep.');

%% Registration test 

[SC_regQ2Q1, warp_Q2Q1] = register_function(grid,Q1,Q2,0.001,1000);
[SC_regQ1Q2, warp_Q1Q2] = register_function(grid,Q2,Q1,0.001,1000);

fprintf('Q2 -> Q1\nOriginal Distance: %0.2f\nRegistered Distance: %0.2f\n\n', ...
    sqrt(sum(diag((Q1 - Q2)*(Q1 - Q2)'))), ...
    sqrt(sum(diag((Q1 - SC_regQ2Q1)*(Q1 - SC_regQ2Q1)'))))

fprintf('Q1 -> Q2\nOriginal Distance: %0.2f\nRegistered Distance: %0.2f\n', ...
    sqrt(sum(diag((Q2 - Q1)*(Q2 - Q1)'))), ...
    sqrt(sum(diag((Q2 - SC_regQ1Q2)*(Q2 - SC_regQ1Q2)'))))

figure(1)
warp_Q1Q2.plot('Q2 -> Q1')
figure(2)
warp_Q2Q1.plot('Q1 -> Q2')