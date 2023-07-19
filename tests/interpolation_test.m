% SCRIPT NAME:
%   interpolation_test
%
% DESCRIPTION:
%   Test the tools used for registration on the spherical manifold 
%
% MATLAB VERSION:
%   R2022b
%

addpath('../')
addpath('../interpolation/')
addpath('../utils/')
addpath('../tangent_basis')

%% Global settings

% resolution settings
ICO_RESOLUTION = 4;

% delta for computing numerical derivatives
delta = 1e-5;

% highest order of spherical harmonics for tangent basis
l = 30; 

%% Setup

% number of vertices on the sphere
P = (4^ICO_RESOLUTION) * 10 + 2;

% generate the coordinates on which to do all our calculations (set up the
% grid) using an icosphere since that is what we'll be working with anyway
ico_mesh = icosphere(ICO_RESOLUTION); % We can load this from a freesurfer mesh instead

% generate a grid object from the mesh
grid = SphericalGrid(ico_mesh, l);

%% Interpolation Test

% generate a test function over the grid
test_fun = @(theta,phi) (sin(theta) .* cos(phi)) * ...
                                (sin(theta) .* cos(phi)).';

grid_data = test_fun(grid.theta, grid.phi);

% generate n points over which to evaluate
n = 200;
theta = linspace(0,pi,n)';
phi = linspace(0,2*pi,n)';

% convert to cartesian coordinates
test_points = sphere_to_cart(theta,phi).';

% setup and run the barycentric interpolation class object
tree = AABBtree(grid);
[Vq,Tq] = tree.get_barycentric_data(test_points,true);
test_values = bary_interp_2D_mex(Vq,Vq,Tq,Tq,grid_data);   
test_values = (test_values + test_values.') / 2;                

% compare interpolated values to the truth
truth = test_fun(theta, phi);
diff = ((truth - test_values).^2);

% display results
fig = figure(1); set(fig,'Units','normalized','OuterPosition',[.05 .3 .85 .45]);
ax1 = subplot(1,3,1);
ax2 = subplot(1,3,2);
ax3 = subplot(1,3,3);
imagesc(truth,'Parent',ax1);
imagesc(test_values,'Parent',ax2);
imagesc(diff,'Parent',ax3);
title(ax1,'Truth')
title(ax2,'Estimated')
title(ax3,'difference')

axis(ax1,'off');
clim(ax1,[-0.3,0.3])
colorbar(ax1)

axis(ax2,'off');
clim(ax2,[-0.3,0.3])
colorbar(ax2)

axis(ax3,'off');
colorbar(ax3)

sgtitle('Interpolation Test')

%% Derivative Test

% generate a test function over the grid
test_fun = @(theta,phi) (sin(theta) .* cos(phi)) * ...
                                (sin(theta) .* cos(phi)).';

grid_data = test_fun(grid.theta, grid.phi);

% the concon object expects a LH and RH warp, so pretend the
% function is in the domain S2 U S2 using a kronecker product 
grid_data = kron(eye(2),grid_data);

% setup and run the barycentric interpolation class object
interp = SphericalInterpolator(grid,grid,delta);

[dxt, dxp] = interp.get_derivative(grid_data);

% the derivative in the direction of the tangent basis
test_dxt = @(theta,phi) (cos(theta) .* cos(phi)) * ...
                                (sin(theta) .* cos(phi)).';
test_dxp = @(theta,phi) (-sin(phi)) * ...
                                (sin(theta) .* cos(phi)).';
test_dyt = @(theta,phi) (sin(theta) .* cos(phi)) * ...
                                (cos(theta) .* cos(phi)).';
test_dyp = @(theta,phi) (sin(theta) .* cos(phi)) * ...
                                (-sin(phi)).';

dxt_mse = mean((dxt(1:P,1:P) - test_dxt(grid.theta, grid.phi)).^2, 'all');
dxp_mse = mean((dxp(1:P,1:P) - test_dxp(grid.theta, grid.phi)).^2, 'all');
dyt_mse = mean((dxt(1:P,1:P)' - test_dyt(grid.theta, grid.phi)).^2, 'all');
dyp_mse = mean((dxp(1:P,1:P)' - test_dyp(grid.theta, grid.phi)).^2, 'all');

fprintf('MSE for spatial derivative estimation: %0.4f, %0.4f, %0.4f, %0.4f,\n', ...
    dxt_mse, dxp_mse, dyt_mse, dyp_mse);

%% Jacobian Test

warp = SphericalWarp(grid,delta);

% test warp and jacobian
warp_theta = grid.theta.^3 ./ pi^2;
warp_phi   = grid.phi;

% manually set the warp vertices
warp.V = sphere_to_cart(warp_theta, warp_phi).';

% trick class into generating the jacobian by composing an identity warp
warp = warp.compose_warp(zeros(P,2));

% the jacobian after mapping to the tangent basis
true_J = abs(((3*grid.theta.^2) ./ pi^2) .* (sin(warp_theta)./(sin(grid.theta)+1e-15)));

% display results
fprintf('MSE for Jacobian estimation: %0.4f\n', mean((warp.J - true_J).^2));
