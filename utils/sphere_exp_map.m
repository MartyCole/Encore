% sphere_exp_map.m - Apply the spherical exponential map in cartestian
% coordinates, from vertices using the tangent vectors gamma
%
% Syntax:  [pts] = sphere_exp_map(vertices,gamma)
%
% Inputs:
%    vertices - Points from which to apply the exponential map
%    gamma - Tangent vector field
%
% Outputs:
%    pts - Points on the unit sphere exp_vertices(gamma)
%

function [pts] = sphere_exp_map(vertices, gamma)
    norm_vec = sqrt(sum(gamma.^2, 2)); 
    idx = norm_vec > eps;
       
    pts = vertices;
    pts(idx,:) = cos(norm_vec(idx)).*pts(idx,:) + ...
                 sin(norm_vec(idx)).*((gamma(idx,:)./norm_vec(idx)));
end