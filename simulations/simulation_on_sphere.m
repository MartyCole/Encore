% SCRIPT NAME:
%   simulation_on_sphere
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
ico_mesh = icosphere(ICO_RESOLUTION); % We can load this from a freesurfer mesh instead

% generate a grid object from the mesh
% NOTE: This takes a little while because I have not yet derived the
% normalisation constants for the tangent basis functions. Currently I am
% numerically estimating the constants over a dense grid.
grid = SphericalGrid(ico_mesh, l);

%% Interpolation Test

% generate a test function over the grid
test_fun = @(theta,phi) (sin(theta) .* cos(phi / 2)) * ...
                                (sin(theta) .* cos(phi / 2)).';

grid_data = test_fun(grid.theta, grid.phi);

% generate n points over which to evaluate
n = 5;
theta = linspace(0,0.4,n)';
phi = linspace(0,0.4,n)';

% convert to cartesian coordinates
test_points = sphere_to_cart(theta,phi).';

% setup and run the barycentric interpolation class object
interpolator = SphericalInterpolator(grid, test_points);
test_values = interpolator.evaluate_2d(grid_data);

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

% setup and run the barycentric interpolation class object
sph_derivative = SphericalDerivative(grid, delta, 'polar');
[dxt, dxp, dyt, dyp] = sph_derivative.get_derivative(grid_data);

test_dxt = @(theta,phi) (cos(theta) .* cos(phi)) * ...
                                (sin(theta) .* cos(phi)).';
test_dxp = @(theta,phi) (sin(theta) .* -sin(phi)) * ...
                                (sin(theta) .* cos(phi)).';
test_dyt = @(theta,phi) (sin(theta) .* cos(phi)) * ...
                                (cos(theta) .* cos(phi)).';
test_dyp = @(theta,phi) (sin(theta) .* cos(phi)) * ...
                                (sin(theta) .* -sin(phi)).';

dxt_mse = mean((dxt' - test_dxt(grid.theta, grid.phi)).^2, 'all');
dxp_mse = mean((dxp' - test_dxp(grid.theta, grid.phi)).^2, 'all');
dyt_mse = mean((dyt' - test_dyt(grid.theta, grid.phi)).^2, 'all');
dyp_mse = mean((dyp' - test_dyp(grid.theta, grid.phi)).^2, 'all');

fprintf('MSE for spatial derivative estimation: %0.4f\n', dxt_mse);

%% Jacobian Test

% test warp and jacobian
warp_theta = grid.theta.^3 ./ pi^2;
warp_phi   = grid.phi;

warp.V = sphere_to_cart(warp_theta, warp_phi).';
warp.J = (3/pi^2) * grid.theta.^2;

% calculate the jacobian of the function using polar coordinates
sph_derivative = SphericalDerivative(grid, delta, 'polar');
est_J = sph_derivative.get_jacobian(warp, false);

% display results
fprintf('MSE for Jacobian estimation: %0.4f\n', mean((warp.J - est_J).^2));

%% Simulate connectivity on the icosphere

N = 500000; % number of fibers to generate
kappa1 = 10; % concentration term (higher means less spread)
kappa2 = 2; % concentration term (higher means less spread)
mu1 = [pi/2,0]; % mean point to create start points around
mu2 = [pi/2,pi]; % mean point to create end points around

mu1_cart = sphere_to_cart(mu1(1),mu1(2));
mu2_cart = sphere_to_cart(mu2(1),mu2(2));

% Generate points from which connections start and end
rng(23072022);
cross_start_points = rand_von_mises(kappa2, 3, mu1_cart, N);
cross_end_points = rand_von_mises(kappa2, 3, mu2_cart, N);

lh_sp = rand_von_mises(kappa2, 3, mu1_cart, N);
lh_ep = rand_von_mises(kappa2, 3, mu1_cart, N);
rh_sp = rand_von_mises(kappa2, 3, mu2_cart, N);
rh_ep = rand_von_mises(kappa2, 3, mu2_cart, N);

mask = (cross_start_points(1,:) < 0) | (cross_end_points(1,:) > 0);
cross_start_points = cross_start_points(:,~mask);
cross_end_points = cross_end_points(:,~mask);

mask = (lh_sp(1,:) < 0) | (lh_ep(1,:) < 0);
lh_sp = lh_sp(:,~mask);
lh_ep = lh_ep(:,~mask);
mask = (rh_sp(1,:) > 0) | (rh_ep(1,:) > 0);
rh_sp = rh_sp(:,~mask);
rh_ep = rh_ep(:,~mask);

start_points = [cross_start_points,lh_sp,rh_sp];
end_points = [cross_end_points,lh_ep,rh_ep];

% display results
fig = figure(2);
title('Subsample of random fiber endpoints')
[X,Y,Z]=sphere(100,100);
surf(0.99*X,0.99*Y,0.99.*Z);

rng(10204382)
idx = randperm(length(start_points), 1000);

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

R = 180/pi; % radius of an idealised spherical brain
FWHM = 30; % bandwidth radius mm
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

% warping the start points
[theta_sp,phi_sp] = cart_to_sphere(start_points');
theta_sp = sqrt(theta_sp)*sqrt(pi);

warped_sp = sphere_to_cart(theta_sp,phi_sp);

% warping the end points          
[theta_ep,phi_ep] = cart_to_sphere(end_points');
theta_ep = sqrt(theta_ep)*sqrt(pi);

warped_ep = sphere_to_cart(theta_ep,phi_ep);

% find closet vertex for each point
interpolator = SphericalInterpolator(grid);
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

%% Interpolation Test 2: Connectivity

% generate the warp (note that this is the inverse to the warp above)
warp = struct('J',[],'V',[]);
warp.V = sphere_to_cart(grid.theta.^2 / pi, grid.phi).';
warp.J = (2*grid.theta) / pi;
warp.J(26)=0;

% setup and run the barycentric interpolation class object
interpolator = SphericalInterpolator(grid, warp.V);
sph_derivative = SphericalDerivative(grid, 1e-5, 'tangent');

% estimate the Jacobian and interpolate the function
warp.J = sph_derivative.get_jacobian(warp, true);
test_SC = interpolator.evaluate_2d(SC) .* (warp.J * warp.J');

figure(5)
imagesc(test_SC)
clim([0,4])
title('Interpolated warped SC')

%% Registration Test

% registration variables
interpolator = SphericalInterpolator(grid);
sph_derivative = SphericalDerivative(grid, 1e-5, 'tangent');
grad_delta = 0.001;
dH = 0;

% cached data
mask = 1-eye(P);
A = grid.A * grid.A';

% initial warp    
warp.J = ones(P,1);
warp.Jmat = zeros(P,4);
warp.Jmat(:,1) = 1;
warp.Jmat(:,4) = 1;
warp.V = grid.V;

foo.J = ones(P,1);
foo.Jmat = zeros(P,4);
foo.Jmat(:,1) = 1;
foo.Jmat(:,4) = 1;
foo.V = grid.V;

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
      
    idx = 1:size(grid.basis,2);%randperm(size(grid.basis,2), 600);
    tmp_dH = -2 * (a * grid.basis(:,idx,1) + b * grid.basis(:,idx,2) + c * grid.laplacian(:,idx));
 
% **** THE FULL COST FUNCTION WITHOUT VECTORISATION ****
%     test = zeros(size(grid.basis,2),1);
%     for i = 1:P
%         for j = 1:P
%             test = test - (2 * (Q1(i,j) - Q2(i,j)) * ( ... 
%                 (dQ2xe1(i,j)*squeeze(grid.basis(i,:,1))' + ...
%                  dQ2xe2(i,j)*squeeze(grid.basis(i,:,2))' + ...
%                  dQ2ye1(i,j)*squeeze(grid.basis(j,:,1))' + ...
%                  dQ2ye2(i,j)*squeeze(grid.basis(j,:,2))' + ...
%                 (0.5*Q2(i,j) * (grid.laplacian(i,:)' + grid.laplacian(j,:)')))) * A(i,j));
%         end
%     end
% **** END FULL COST FUNCTION WITHOUT VECTORISATION ****

    if true%norm(tmp_dH,'fro') > norm(dH,'fro')
        dH = tmp_dH;

        % calculate displacement in each basis direction
        test_gamma = squeeze(sum(-dH.*grid.basis(:,idx,:),2));
        
        % compose the new warp with all previous warps
        [warp,foo] = sph_derivative.compose_warp(grad_delta .* test_gamma, foo);
        foo.J = sph_derivative.get_jacobian(foo, true);  

        warped = true;
    else
        fprintf('Converged (increased norm) %d: %0.6f\n', iter, cost)            
        break   
    end    

    % evaluate function after warping
    interpolator.update_query_points(foo.V);
    moving_img = interpolator.evaluate_2d(Q2) .* (sqrt(foo.J) * sqrt(foo.J).');
    moving_img = (moving_img + moving_img.') / 2;
    moving_img = moving_img / sqrt(sum(moving_img(:).^2 .* A(:)));
    
    % evaluate the new cost
    FmM = Q1 - moving_img; 
    cost = sum(FmM(:).^2 .* A(:));
   
    if ((last_cost - cost) < -1e-2)       
       fprintf('Converged (increased cost) %d: %0.6f\n', iter, cost)            
       break
    end

    last_cost = cost; 

    if (mod(iter,10) == 0)
        pause(0.5);
        fprintf('Iteration %d cost: %0.6f, %0.4f\n', iter, cost, sum(moving_img(:).^2 .* A(:)))
        figure(6)
        foo.J = abs((foo.Jmat(:,1).*foo.Jmat(:,4)) - (foo.Jmat(:,2).*foo.Jmat(:,3)));
        trisurf(grid.T,foo.V(:,1),foo.V(:,2),foo.V(:,3),foo.J)
        axis off   
        axis equal
        caxis([0,2]);
        pause(0.5);
    end
end   

%% Display the final warp

figure(6)
foo.J = abs((foo.Jmat(:,1).*foo.Jmat(:,4)) - (foo.Jmat(:,2).*foo.Jmat(:,3)));
trisurf(grid.T,foo.V(:,1),foo.V(:,2),foo.V(:,3),foo.J)
axis off
axis equal