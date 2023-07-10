% SCRIPT NAME:
%   sphere_registration_test
%
% DESCRIPTION:
%   Test the tools used for registration on the spherical manifold (just
%   the one hemisphere) and perform a test registration in both directions,
%   Q1 -> Q2 and Q2 >- Q1. Tools tested are Barycentric Interpolation, and
%   Jacobian estimation.
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

% maximum number of iterations for registration
REG_ITERS = 40000;

%% Setup

% number of vertices on the sphere
P = (4^ICO_RESOLUTION) * 10 + 2;

% generate the coordinates on which to do all our calculations (set up the
% grid) using an icosphere since that is what we'll be working with anyway
ico_mesh = icosphere(ICO_RESOLUTION); 

% generate a grid object from the mesh
grid = SphericalGrid(ico_mesh, l);
sph_derivative = SphericalDerivative(grid, 1e-5, 'tangent');

%% Simulate diffeomorphism

% Set seed for reproduciblity
rng(23711228);

% initialise an identity warp
test_warp = SphericalWarp(grid);

% simulate exponential decay
decay = exp([linspace(0,-8,grid.num_basis/2),linspace(0,-8,grid.num_basis/2)]);

% generate coefficients for all directions
coeff = (2 .* rand(1,grid.num_basis) - 1) .* decay;
coeff(randperm(grid.num_basis, floor(grid.num_basis .* 0.75))) = 0;

% generate the displacement vector
gamma = squeeze(sum(-coeff .* grid.basis,2));

% compose the warp multiple times
for i = 1:8000
    test_warp = sph_derivative.compose_warp(0.00001 .* gamma, test_warp);
end

% generate the Jacobian
test_warp.J = sph_derivative.get_jacobian(test_warp, true);  

figure(1)
title('Simulated Diffeomorphism')
trisurf(grid.T,test_warp.V(:,1),test_warp.V(:,2),test_warp.V(:,3),test_warp.J)
axis off
axis equal

%% Simulate connectivity on the icosphere

% Set seed for reproduciblity
rng(23072022);

N = 500000;       % number of fibers to generate
kappa1 = 10;      % concentration term (higher means less spread)
kappa2 = 2;       % concentration term (higher means less spread)
mu1 = [pi/2,0];   % mean point to create start points around
mu2 = [pi/2,pi];  % mean point to create end points around

mu1_cart = sphere_to_cart(mu1(1),mu1(2));
mu2_cart = sphere_to_cart(mu2(1),mu2(2));

% Generate points from which connections start and end
lr_sp = rand_von_mises(kappa2, 3, mu1_cart, N);
lr_ep = rand_von_mises(kappa2, 3, mu2_cart, N);
lh_sp = rand_von_mises(kappa2, 3, mu1_cart, N);
lh_ep = rand_von_mises(kappa2, 3, mu1_cart, N);
rh_sp = rand_von_mises(kappa2, 3, mu2_cart, N);
rh_ep = rand_von_mises(kappa2, 3, mu2_cart, N);

% Mask LH points to be only on LH
mask = (lh_sp(1,:) < 0) | (lh_ep(1,:) < 0);
lh_sp = lh_sp(:,~mask);
lh_ep = lh_ep(:,~mask);

% Mask RH points to be only on RH
mask = (rh_sp(1,:) > 0) | (rh_ep(1,:) > 0);
rh_sp = rh_sp(:,~mask);
rh_ep = rh_ep(:,~mask);

% Concatenate all points
start_points = [lr_sp,lh_sp,rh_sp];
end_points = [lr_ep,lh_ep,rh_ep];

% display results
fig = figure(2);
title('Subsample of random fiber endpoints')
[X,Y,Z]=sphere(100,100);
surf(0.99*X,0.99*Y,0.99.*Z);

% get a subsample of points to keep things neat
idx = randperm(length(start_points), 5000);

hold on;
shading flat;
colormap([.2 0 .8]);
scatter3(start_points(1,idx),start_points(2,idx),start_points(3,idx),...
    'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.5,'MarkerFaceColor',[1,0,0]);
scatter3(end_points(1,idx),end_points(2,idx),end_points(3,idx),...
    'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.5,'MarkerFaceColor',[0,1,0]);
axis equal;
axis off;
hold off;
view(0,0);

%% Calculate Kernel for KDE smoothing

R = 180/pi;    % radius of an idealised spherical brain
FWHM = 15;     % bandwidth radius mm
epsilon = 0.1; % truncation distance

% get geodesic distance between each pair of vertices on the grid
dst = zeros(P,P);

for i = 1:P
    for j = (i+1):P
        dst(i,j) = acos(max(min(dot(grid.V(i,:), grid.V(j,:)),1),-1))*R;
    end
end

dst = dst + dst';

% calculate the bandwidth and truncation distance
sigma = (FWHM/2) / sqrt(8*log(2));
max_dist = (FWHM/2) * sqrt((-log2(epsilon)));

% calculate the kernel for the grid
kernel = exp(-(dst.^2 / (2 * (sigma^2))));
kernel(dst > max_dist) = 0;

% normalise the kernel so that counts stay the same
kernel = kernel ./ sum(kernel,2);

%% Calculate KDE for the original endpoints

% find closet vertex for each point
interpolator = SphericalInterpolator(grid);
[D_sp,T_sp] = interpolator.get_barycentric_data(start_points');
[D_ep,T_ep] = interpolator.get_barycentric_data(end_points');

[~,idx_a] = min(D_sp,[],2,'linear');
[~,idx_b] = min(D_ep,[],2,'linear');

% create adjacency matrix from endpoints
node_a = T_sp(idx_a);
node_b = T_ep(idx_b);
A = full(sparse(node_a, node_b, 1, P, P));
A = A + A';

% generate smooth SC
SC = (kernel*A*kernel');

% display results
figure(3)
imagesc(SC)
clim([0,4])
title('Simulated SC')

%% Calculate KDE for the warped endpoints

interpolator = SphericalInterpolator(grid);
[D_sp,T_sp] = interpolator.get_barycentric_data(start_points');
[D_ep,T_ep] = interpolator.get_barycentric_data(end_points');

wx_sp = dot(D_sp, reshape(test_warp.V(T_sp,1),length(start_points),3),2);
wy_sp = dot(D_sp, reshape(test_warp.V(T_sp,2),length(start_points),3),2);
wz_sp = dot(D_sp, reshape(test_warp.V(T_sp,3),length(start_points),3),2);            

warped_sp = normr([wx_sp,wy_sp,wz_sp]).';

wx_ep = dot(D_ep, reshape(test_warp.V(T_ep,1),length(start_points),3),2);
wy_ep = dot(D_ep, reshape(test_warp.V(T_ep,2),length(start_points),3),2);
wz_ep = dot(D_ep, reshape(test_warp.V(T_ep,3),length(start_points),3),2);            

warped_ep = normr([wx_ep,wy_ep,wz_ep]).';

% find closet vertex for each point
[D_sp,T_sp] = interpolator.get_barycentric_data(warped_sp');
[D_ep,T_ep] = interpolator.get_barycentric_data(warped_ep');

[~,idx_a] = min(D_sp,[],2,'linear');
[~,idx_b] = min(D_ep,[],2,'linear');

% create adjacency matrix from endpoints
node_a = T_sp(idx_a);
node_b = T_ep(idx_b);
A = full(sparse(node_a, node_b, 1, P, P));
A = A + A';

% generate warped SC
warped_SC = (kernel*A*kernel');

% display results
figure(4)
imagesc(warped_SC)
clim([0,4])
title('True warped SC')

%% Registration Test

% registration variables
interpolator = SphericalInterpolator(grid);
sph_derivative = SphericalDerivative(grid, 1e-5, 'tangent');

grad_delta = 0.001;

% cached data
mask = 1-eye(P);
A = grid.A * grid.A';

% initial warp    
warp = SphericalWarp(grid);

% initial functions
Q1 = sqrt(SC / sum(SC(:) .* A(:)));
Q2 = sqrt(warped_SC / sum(warped_SC(:) .* A(:)));

% initial cost
moving_img = Q2;
FmM = (Q1 - moving_img);
last_cost = sum(FmM(:).^2 .* A(:));

fprintf('Initial cost for Sub: %0.6f\n', last_cost)
init_cost = last_cost;

% start registering
for iter = 1:REG_ITERS    
    warped = false;

    % calculate the derivative
    [dQ2xe1, dQ2xe2, dQ2ye1, dQ2ye2] = sph_derivative.get_derivative(moving_img);

    % evaluate derivative of the cost function
    FmM = FmM .* A;

    a = sum(FmM .* ((dQ2xe1' + dQ2ye1)), 1);
    b = sum(FmM .* ((dQ2xe2' + dQ2ye2)), 1);
    c = sum(FmM .* moving_img, 1);
      
    dH = -2 * ((a * grid.basis(:,:,1)) + (b * grid.basis(:,:,2)) + (c * grid.laplacian));

    if norm(dH,'fro') < 1e-6
        fprintf('Converged (small norm) %d: %0.6f\n', iter, cost)            
        break   
    end    

    % calculate displacement in each basis direction
    gamma = squeeze(sum(-dH .* grid.basis,2));            
    step_size = min(0.02 / max(vecnorm(gamma')), grad_delta);

    % compose the new warp with all previous warps
    warp = sph_derivative.compose_warp(step_size .* gamma, warp);
    warp.J = sph_derivative.get_jacobian(warp,true);
    
    % evaluate function after warping
    interpolator.update_query_points(warp.V);
    moving_img = interpolator.evaluate_2d(Q2) .* (sqrt(warp.J) * sqrt(warp.J).');
    moving_img = (moving_img + moving_img.') / 2;
    moving_img = moving_img / sqrt(sum(moving_img(:).^2 .* A(:)));
    
    % evaluate the new cost
    FmM = Q1 - moving_img; 
    cost = sum(FmM(:).^2 .* A(:));
   
    if ((last_cost - cost) < -1e-3)       
       fprintf('Converged (increased cost) %d: %0.6f\n', iter, cost)            
       break
    end

    last_cost = cost; 

    if (mod(iter,10) == 0)       
        fprintf('Iteration %d cost: %0.6f, %0.4f\n', iter, cost, sum(moving_img(:).^2 .* A(:)));       
    end
end   

%% Display the final warp

figure(7)
title('Estimated Warp')
trisurf(grid.T,warp.V(:,1),warp.V(:,2),warp.V(:,3),warp.J)
axis off
axis equal
