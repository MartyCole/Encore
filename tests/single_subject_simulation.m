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
lh_warp = random_diffeomorphism(grid,8,0.75,8000);
rh_warp = random_diffeomorphism(grid,8,0.75,8000);

% generate the connectomes
[orig_sp,orig_ep] = random_connectome(N,kappa1,mu1,kappa2,mu2);

concon = Concon(grid,grid,orig_sp.',orig_ep.',10,0.1,1e-5);
F1 = concon.evaluate();
F2 = concon.evaluate(lh_warp,rh_warp);

%% Registration test 
 
encore = Encore(grid,grid,33,0.01,1000,0); 
[reg_lh,reg_rh,~] = encore.register(F1,F2);

figure(1); lh_warp.plot('LH warp'); 
figure(2); reg_lh.plot('LH reg'); 
figure(3); rh_warp.plot('RH warp'); 
figure(4); reg_rh.plot('RH reg');
