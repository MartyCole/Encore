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
F1 = spherical_kde(grid, kernel, orig_sp.', orig_ep.');
F2 = spherical_kde(grid, kernel, warp_sp.', warp_ep.');

Q1 = [F1,zeros(size(F1));zeros(size(F1)),F2];
Q2 = [F2,zeros(size(F1));zeros(size(F1)),F1];

%% Registration test 

encore = Encore(grid,grid,30,0.001,1000,0);

[SC_reg,lh_warp,rh_warp,~] = encore.register(Q1,Q2);

fprintf('Q2 -> Q1\nOriginal Distance: %0.2f\nRegistered Distance: %0.2f\n\n', ...
     sqrt(sum(diag((Q1 - Q2)*(Q1 - Q2)'))), ...
     sqrt(sum(diag((Q1 - SC_reg)*(Q1 - SC_reg)'))))
 
figure(1)
lh_warp.plot('Q2 -> Q1')
figure(2)
rh_warp.plot('Q1 -> Q2')
figure(3)
imagesc(Q1 + [zeros(size(F1)), SC_reg(1:P,1:P); SC_reg((P+1):end,(P+1):end), zeros(size(F1));])
