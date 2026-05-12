% sphere_log_map.m - Apply the spherical logarithm map in cartestian
% coordinates.
%
% Syntax:  [pts] = sphere_log_map(mu,q)
%
% Inputs:
%    mu - Points from which to apply the log map
%    q - Target points to lift from the unit sphere
%
% Outputs:
%    pts - Tangent vectors log_mu(q)
%

function [gamma] = sphere_log_map(mu, q)
    theta = acos(max(-1, min(1, dot(q,mu,2))));
    idx = theta > eps;                
            
    gamma = zeros(size(mu));
    gamma(idx,:) = (theta(idx) ./ sin(theta(idx))) .* (q(idx,:) - (dot(q(idx,:),mu(idx,:),2) .* mu(idx,:)));   
end