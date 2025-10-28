% MATLAB R2018a
%
% FUNCTION NAME:
%   simulate_connectivity1D
%
% DESCRIPTION:
%   Simulate a pxp connectivity matrix for the 1 dimensional domain [0 1]. 
%
% INPUT:
%   N - Integer with the number of connections to generate
%   p - Integer for the number of regions of interest
%   mu1 - Scalar indicating the mean value for connection start points
%   mu2 - Scalar indicating the mean value for connection end points
%   sd - Scalar for the standard deviation of connection points
%   circluar - Boolean indicating if the domain is on a line or a circle
%
% OUTPUT:
%   C - pxp matrix of connection counts within each region of interest
%   start_points - Nx1 vector of values between [0,1] of the start points
%   end_points - Nx1 vector of values between [0,1] of the end points
%   Side effects: none
%
% ASSUMPTIONS AND LIMITATIONS:
%   Points are distributed normally on the domain

function [C, start_points, end_points] = simulate_connectivity1D(N, p, mu1, mu2, sd, circular)

% generate points from which connections start
start_points = [normrnd(mu1, sd, N, 1)]; 
                %normrnd(mu2, sd/2, floor(sqrt(N)), 1)
                %normrnd(mu1, sd/2, floor(sqrt(N)), 1)];

% generate points where connections end
end_points = [normrnd(mu2, sd, N, 1)]; 
              %normrnd(mu2, sd/2, floor(sqrt(N)), 1)
              %normrnd(mu1, sd/2, floor(sqrt(N)), 1)]; 

if circular == true
    start_points = mod(start_points, 1);
    end_points = mod(end_points, 1);
else
    % truncate connections to withiin [0, 1]
    in_domain_1 = (start_points > 0) & (start_points < 1);
    in_domain_2 = (end_points > 0) & (end_points < 1);

    % find pairs of points that are fully within [0, 1]
    in_domain = in_domain_1 & in_domain_2;

    start_points = start_points(in_domain);
    end_points = end_points(in_domain);    
end

sp = [start_points;end_points];
ep = [end_points;start_points];

C = accumarray([round((sp(:) * (p-1)) + 1), round((ep(:) * (p-1)) + 1)], 1, [p,p]);

end