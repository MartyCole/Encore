% SCRIPT NAME:
%   SphericalWarp
%
% DESCRIPTION:
%   Class to represent a diffeomorphism on the sphere.
%
% MATLAB VERSION:
%   R2022b
%
classdef SphericalWarp
    properties
        V
        J
        Jdet
    end
    methods
        function obj = SphericalWarp(grid)
            P = size(grid.V,1);

            % The warped vertices
            obj.V = grid.V;

            % Determinate of the Jacobian Matrix
            obj.Jdet = ones(P,1);
            
            % Jacobian matrix at each vertex
            obj.J = zeros(P,4);
            obj.J(:,1) = 1;
            obj.J(:,4) = 1;           
        end
    end
end