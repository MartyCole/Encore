% sphere_log_map.m - Apply the spherical logarithm map in cartestian
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

function [gamma] = sphere_log_map(mu, q)
    tmp_theta = acos(dot(q,mu,2));
    idx = tmp_theta > eps;                
            
    gamma = zeros(size(mu));
    gamma(idx,:) = (tmp_theta(idx) ./ sin(tmp_theta(idx))) .* (q(idx,:) - (dot(q(idx,:),mu(idx,:),2) .* mu(idx,:)));   
end