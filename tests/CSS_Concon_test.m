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

rng(23711228); % set seed for reproduciblity
l = 31;        % order of SH for tangent basis

%% Simulate connectivity on the icosphere

% load real data example
xing = readmatrix("../data/100206_xing_sphere_avg_coords.tsv",'FileType','text');

mesh = read_vtk("../data/lh_sphere_avg_0.94.vtk");
mesh.V = normr(mesh.vtx.');
mesh.T = mesh.tri;
lh_grid = SphericalGrid(mesh, l);

mesh = read_vtk("../data/rh_sphere_avg_0.94.vtk");
mesh.V = normr(mesh.vtx.');
mesh.T = mesh.tri;
rh_grid = SphericalGrid(mesh, l);

hemi1 = 1-xing(:,2);
hemi2 = 1-xing(:,7);

sp = xing(:,3:5)';
ep = xing(:,8:10)';

%% Calculate the heat kernels

% calculate smoothed connectomes
hcp_kernel = spherical_heat_kernel(lh_grid,rh_grid,0.005,l);

moyer_sc = load("../data/smoothed_sc_avg_0.005_0.94.mat", 'sc');
moyer_sc = moyer_sc.sc + moyer_sc.sc.';
moyer_sc = moyer_sc ./ sum(moyer_sc(:));

css_concon = Concon(lh_grid, rh_grid, sp.',ep.',hemi1,hemi2);
css_sc = css_concon.evaluate(hcp_kernel, true);

fig = figure(1);
t = tiledlayout(1,2);
t.Padding = 'compact';
t.TileSpacing = 'compact';

nexttile
imagesc(log(1e6*moyer_sc + 1));
clim([0,6])
xticklabels([]);
yticklabels([]);
title('Moyer SC')
nexttile
imagesc(log(1e6*css_sc + 1));
clim([0,6])
xticklabels([]);
yticklabels([]);
colorbar()
title('CSS SC')

%% Downsample Test

ico_mesh = icosphere(4);                 % generate a triagulated mesh for the grid
grid = SphericalGrid(ico_mesh, l);       % generate a grid object from the mesh

lo_ico_mesh = icosphere(3);              % generate a triagulated mesh for the grid
lo_grid = SphericalGrid(lo_ico_mesh, l); % generate a grid object from the mesh

% calculate smoothed connectomes
kernel = spherical_heat_kernel(grid,grid,0.005,l);
lo_kernel = spherical_heat_kernel(lo_grid,lo_grid,0.005,l);

F1 = Concon(grid, grid, sp.',ep.',hemi1,hemi2);
F2 = copy(F1);
F2.resample(lo_grid,lo_grid);

fig = figure(2);
t = tiledlayout(1,2);
t.Padding = 'compact';
t.TileSpacing = 'compact';

nexttile
imagesc(log(1e6*F1.evaluate(kernel,true) + 1));
clim([0,6])
xticklabels([]);
yticklabels([]);
title('CSS SC')
nexttile
imagesc(log(1e5*F2.evaluate(lo_kernel,true) + 1));
clim([0,6])
xticklabels([]);
yticklabels([]);
colorbar()
title('Resampled CSS SC')

