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
        T
    end
    methods
        function obj = SphericalWarp(grid)
            P = size(grid.V,1);

            % The warped vertices
            obj.V = grid.V;

            % Determinate of the Jacobian Matrix
            obj.J = ones(P,1);       

            % The triangulation
            obj.T = grid.T;
        end

        function plot(obj, fig_title)            
            title(fig_title)
            trisurf(obj.T, obj.V(:,1), obj.V(:,2), obj.V(:,3), obj.J)
            axis off
            axis equal
        end
    end
end