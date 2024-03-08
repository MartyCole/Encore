% MATLAB R2018a
%
% FUNCTION NAME:
%   warp_connectivity1D
%
% DESCRIPTION:
%   Uses warp on start_points and end_points to generate a pxp connectivity 
%   matrix. Also returns the warped start and end points. 
%
% INPUT:
%   N - Integer with the number of connections to generate
%   p - Integer for the number of regions of interest
%   mu1 - Scalar indicating the mean value for connection start points
%   mu2 - Scalar indicating the mean value for connection end points
%   sd - Scalar for the standard deviation of connection points
%
% OUTPUT:
%   C - pxp matrix of connection counts within each region of interest
%   start_points - Nx1 vector of values between [0,1] of the start points
%   end_points - Nx1 vector of values between [0,1] of the end points
%   Side effects: none
%
% ASSUMPTIONS AND LIMITATIONS:
%   Points are distributed normally on the domain

function [C, new_start_point, new_end_point] = warp_connectivity1D(start_point, end_point, warp, p)

warp = invert_gamma(warp')';
warp = (warp * (p-1)) + 1;

start_point = (start_point * (length(warp)-1)) + 1;
end_point = (end_point * (length(warp)-1)) + 1;

idx1 = floor(start_point);
idx2 = floor(end_point);
area1 = mod(start_point, 1);
area2 = mod(end_point, 1);

new_start_point = warp(idx1)+area1.*(warp(idx1+1)-warp(idx1));
new_end_point = warp(idx2)+area2.*(warp(idx2+1)-warp(idx2));

% generate the connectivity matrix
sp = round([new_start_point;new_end_point]);
ep = round([new_end_point;new_start_point]);

C = accumarray([sp(:), ep(:)], 1, [p,p]);

end