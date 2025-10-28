% SCRIPT NAME:
%   SphericalGrid
%
% DESCRIPTION:
%   Class to represent a grid over which to registration. The input is a 
%   triangular mesh of a sphere.
%
% MATLAB VERSION:
%   R2022b
%
classdef SphericalGrid
    properties
        V
        T
        theta
        phi
        A 
        basis
        laplacian
        e1
        e2    
        num_basis
    end
    methods
        function obj = SphericalGrid(mesh, l)
            % normalise mesh to unit sphere
            obj.V = normr(mesh.V);
            obj.T = mesh.T;
            
            % convert cartesian to spherical coordinates
            [obj.theta, obj.phi] = cart_to_sphere(obj.V);

            % calculate voronoi areas of the grid
            obj.A = calc_voronoi_area(obj.V, obj.T);
            
            % generate the basis system at points on the icosphere           
            [obj.basis, obj.laplacian, ~, ~, ~] = tangent_basis(l, obj.theta, obj.phi, obj.A);
            obj.basis = squeeze(obj.basis);
            obj.laplacian = squeeze(obj.laplacian);
           
            % save the number of basis functions
            obj.num_basis = size(obj.basis,2);

            % basis system for cartesian tangent space
            obj.e1 = [cos(obj.theta).*cos(obj.phi), ...
                      cos(obj.theta).*sin(obj.phi), ...
                     -sin(obj.theta)];
    
            obj.e2 = [-sin(obj.phi), ...
                       cos(obj.phi), ...
                       zeros(size(obj.theta,1),1)];
        end
    end
end