function kernel = spherical_kernel(grid,FWHM,epsilon)

% radius of an idealised spherical brain
R = 180/pi;    

% get geodesic distance between each pair of vertices on the grid
dst = acos(max(min(grid.V * grid.V.',1),-1))*R;

% calculate the bandwidth and truncation distance
sigma = (FWHM/2) / sqrt(8*log(2));
max_dist = (FWHM/2) * sqrt((-log2(epsilon)));

% calculate the kernel for the grid
kernel = exp(-(dst.^2 / (2 * (sigma^2))));
kernel(dst > max_dist) = 0;

% normalise the kernel so that counts stay the same
kernel = kernel ./ sum(kernel,2);
kernel = kron(eye(2),kernel);

end

