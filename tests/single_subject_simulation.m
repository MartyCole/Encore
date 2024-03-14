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
[orig_sp1,orig_ep1] = random_connectome(floor(N*0.75),1,[pi/2,0],2,[pi/2,pi]);
[orig_sp2,orig_ep2] = random_connectome(floor(N*0.25),3,[pi/4,pi/3],3,[pi/3,pi/4]);

orig_sp = [orig_sp1 orig_sp2];
orig_ep = [orig_ep1 orig_ep2];

clear("orig_sp1","orig_ep1","orig_sp2","orig_ep2")

hemi1 = randsample([0 1], size(orig_ep,2), true);
hemi2 = randsample([0 1], size(orig_ep,2), true);
hemi2(1:(0.8*size(orig_ep,2))) = hemi1(1:(0.8*size(orig_ep,2)));

% calculate smoothed connectomes
kernel = spherical_kernel(grid,grid,5,0.1);

F1 = Concon(grid,grid,orig_sp.',orig_ep.',hemi1,hemi2);
F2 = copy(F1);
F2.warp_connectome(lh_warp,rh_warp);

%% Registration test 

encore = Encore(grid,grid,30,0.005,1000,0);

[Freg,est_lh_warp,est_rh_warp,cost] = encore.register(F2,F1,kernel,'verbose',1);

%% Results

orig_F1 = F1.evaluate(kernel, true);
orig_F2 = F2.evaluate(kernel, true);
reg_F1toF2 = Freg.evaluate(kernel, true);

fprintf('Q2 -> Q1\nOriginal Distance: %0.2f\nRegistered Distance: %0.2f\n\n', ...
     sqrt(sum(diag((orig_F1 - orig_F2)*(orig_F1 - orig_F2)'))), ...
     sqrt(sum(diag((reg_F1toF2 - orig_F2)*(reg_F1toF2 - orig_F2)'))))
 
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

