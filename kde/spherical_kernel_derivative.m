function kernel = spherical_kernel_derivative(lh_grid,rh_grid,FWHM,epsilon)

% radius of an idealised spherical brain
R = 180/pi;    

% calculate the bandwidth and truncation distance
sigma = (FWHM/2) / sqrt(8*log(2));
max_dist = (FWHM/2) * sqrt((-log2(epsilon)));

grid{1} = lh_grid;
grid{2} = rh_grid;

kernel = cell(2,1);

for i = 1:2
    % get geodesic distance between each pair of vertices on the grid
    dst = acos(max(min(grid{i}.V * grid{i}.V.',1),-1))*R;

    % calculate the kernel for the grid    
    kernel{i} = (dst.^2/(2*pi*sigma^6)) .* exp(-(dst.^2 / ((sigma^2))));   
    kernel{i}(dst > max_dist) = 0;

    % normalise the kernel so that counts stay the same
    kernel{i} = kernel{i} ./ sum(kernel{i},1);
end

lh_N = size(lh_grid.V,1);
rh_N = size(rh_grid.V,1);

kernel = [kernel{1} zeros(lh_N,rh_N); 
          zeros(rh_N,lh_N) kernel{2}];

end

