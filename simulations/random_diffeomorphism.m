function warp = random_diffeomorphism(grid,decayrate,density,iters)

% initialise an identity warp
warp = SphericalWarp(grid);

% simulate exponential decay
decay = exp([linspace(0,-decayrate,grid.num_basis/2), ...
             linspace(0,-decayrate,grid.num_basis/2)]);

% generate coefficients for all directions
coeff = (2 .* rand(1,grid.num_basis) - 1) .* decay;
coeff(randperm(grid.num_basis, floor(grid.num_basis .* density))) = 0;

% generate the displacement vector
gamma = squeeze(sum(-coeff .* grid.basis, 2));

% compose the warp multiple times
sph_derivative = SphericalDerivative(grid, 1e-5, 'tangent');

for i = 1:iters
    warp = sph_derivative.compose_warp(0.00001 .* gamma, warp);
end

% generate the Jacobian
warp.J = sph_derivative.get_jacobian(warp, true);  

end

