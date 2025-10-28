% MATLAB R2018a
%
% FUNCTION NAME:
%   random_warp1D
%
% DESCRIPTION:
%   Find the point and index where a ray intersects one of the triangles
%   composed of the vertices in the vectors [vertex1, vertex2, vertex3]. 
%
% INPUT:
%   ray - 1x3 normalised vector in the direction of the ray
%   vertex1 - Nx3 matrix of the first vertices of the N triangles
%   vertex2 - Nx3 matrix of the second vertices of the N triangles
%   vertex3 - Nx3 matrix of the third vertices of the N triangles
%
% OUTPUT:
%   p - 1x3 vector of the intersection point
%   idx - Nx1 boolean vector indicating the index of the intersecting
%         triangle
%   Side effects: none
%
% ASSUMPTIONS AND LIMITATIONS:
%   The ray is assumed to originate at [0, 0, 0]

function warp = random_noisy_warp1D(N, T, sigma)

warp = zeros(N,T);
sample_points = ceil(sqrt(T));
warp_points = rand(N,sample_points).^sigma;

for k = 1:N
    tmp_warp = cumtrapz(warp_points(k,:));
    tmp_warp = tmp_warp / tmp_warp(end);
    
    warp(k,:) = interp1(linspace(1,T,sample_points), tmp_warp, 1:T, 'pchip');
end

end