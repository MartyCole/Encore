% SCRIPT NAME:
%   real_data_simulation
%
% DESCRIPTION:
%   Loads real data, warps the connectome and recovers the original warp
%   through registration. Demonstrates Q1->Q2 and Q2->Q1 have same accuracy
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

%% Load connectivity on the icosphere

% load real data example (endpoints on a sphere)
xing = readmatrix("../data/100206_xing_sphere_avg_coords.tsv",'FileType','text');

start_hemi = 1-xing(:,2); % original Concon had the left hemisphere coded as 1
end_hemi = 1-xing(:,7);   % original Concon had the left hemisphere coded as 1

start_pts = xing(:,3:5);
end_points = xing(:,8:10);

% generate a random warp
lh_warp = random_diffeomorphism(grid,8,0.75,8000);
rh_warp = random_diffeomorphism(grid,8,0.75,8000);

% calculate heat kernel for smoothing connectomes
kernel = spherical_heat_kernel(grid,grid,0.005,l);

%% Registration test 
F1 = Concon(grid,grid,start_pts,end_points,start_hemi,end_hemi);
F2 = copy(F1);
F2.warp_connectome(lh_warp,rh_warp);
F2.resample(grid,grid); % resets warp to identity for registration

% the registration algorithm class object (note that smaller l encourages
% smoother diffeomorphisms due to the simpler number of directions)
encore = Encore(ico_mesh,ico_mesh,l,0.05,1000,1e-5);

% register in both directions
[F1_to_F2,~,~,~] = encore.register(F2,F1,kernel,'verbose',1);
[F2_to_F1,~,~,~] = encore.register(F1,F2,kernel,'verbose',2);

%% Display the original warps for comparison
 
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

%% Show the results

orig_F1 = F1.evaluate(kernel, true);
orig_F2 = F2.evaluate(kernel, true);
reg_F1toF2 = F1_to_F2.evaluate(kernel, true);
reg_F2toF1 = F2_to_F1.evaluate(kernel, true);

fprintf('Q2 -> Q1\nOriginal Distance: %0.2f\nF1->F2 Distance: %0.2f\nF2->F1 Distance: %0.2f\n\n', ...
     acos(sum(sqrt(orig_F1).*sqrt(orig_F2),'all')), ...
     acos(sum(sqrt(reg_F1toF2).*sqrt(orig_F2),'all')), ...
     acos(sum(sqrt(reg_F2toF1).*sqrt(orig_F1),'all')))

fig = figure(4);
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
