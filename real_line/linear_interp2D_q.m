% MATLAB R2018a
%
% FUNCTION NAME:
%   linear_interp2D
%
% DESCRIPTION:
%   Performs bilinear interpolation of a connectivity function on the
%   domain of [0, 1]. Note that this is not the same as warping the
%   connections themselves.
%
% INPUT:
%
% OUTPUT:
%   C - pxp matrix of connection counts within each region of interest
%   start_points - Nx1 vector of values between [0,1] of the start points
%   end_points - Nx1 vector of values between [0,1] of the end points
%   Side effects: none
%
% ASSUMPTIONS AND LIMITATIONS:
%   Points are distributed normally on the domain

function new_C = linear_interp2D_q(C, warp)

p = size(C, 1);
stepsize = 1 / (p-1);

grad = sqrt(gradient(warp', stepsize)') * sqrt(gradient(warp', stepsize)')';

x = repmat(linspace(0,1,p),p,1);
y = repmat(linspace(0,1,p)',1,p);

new_x = repmat(warp',p,1);
new_y = repmat(warp,1,p);
new_C = interp2(x,y,C,new_x,new_y,'linear') .* grad;

end