classdef QuarticKernel < Kernel
    % QuarticKernel
    %
    % Implements the quartic (biweight) kernel on the unit sphere S^2
    %
    % This class provides:
    %   - fast evaluation of the quartic kernel matrix K for a given bandwidth σ
    %   - directional derivatives of the kernel (x, y, z components)
    %
    % Internally, the kernel is computed using:
    %       K = (15/16σ) {1 - (X/σ)^2}^2 for X ≤ σ
    %
    % This class subclasses Kernel and implements compute_sigma.

    properties (Access = private)
        lh_X
        rh_X        
    end

    methods
        function obj = QuarticKernel(lh_grid, rh_grid, gpu)
            % Construct a quartic kernel object.
            %
            % Inputs:
            %   lh_grid, rh_grid : SphericalGrid objects for LH and RH          
            %   gpu              : (true) whether to store kernel matrices on GPU                   

            if nargin < 3
                gpu = true;
            end

            obj.use_GPU = gpu;

            % save the grid the kernel is on
            obj.lh_grid = lh_grid;
            obj.rh_grid = rh_grid;           

            % cache the distance between vertices on the grid
            obj.lh_X = obj.calc_distance(lh_grid.V * lh_grid.V.');
            obj.rh_X = obj.calc_distance(rh_grid.V * rh_grid.V.');
        end      

        function [K, dK] = compute_sigma(obj, sigma)
            % Compute the quartic kernel matrix K for bandwidth sigma.        
            %
            % Outputs:
            %   K  : kernel matrix
            %   dK : struct of directional derivatives (if requested)    

            % compute the kernel with the given sigma
            lh_K = (15/16) * (1/sigma) * (1 - (obj.lh_X / sigma)^2)^2;
            rh_K = (15/16) * (1/sigma) * (1 - (obj.rh_X / sigma)^2)^2;

            % only support within sigma
            lh_K(obj.lh_X > sigma) = 0;
            rh_K(obj.rh_X > sigma) = 0;
            
            K = blkdiag(lh_K, rh_K);            

            % row normalise the kernel (diffusion smoothing)
            K_norm = sum(K,2);   
            K = K ./ K_norm;

            % move to the GPU?
            K = obj.maybe_GPU(K);

            % compute the directional derivative if requested
            if nargout == 2                
                % compute the derivative of the kernel with the given sigma
                lh_dK = -(15/4) * (1/sigma^4) * (1 - (obj.lh_X^2 / sigma^2));
                rh_dK = -(15/4) * (1/sigma^4) * (1 - (obj.rh_X^2 / sigma^2));

                % only support within sigma
                lh_dK(obj.lh_X > sigma) = 0;
                rh_dK(obj.rh_X > sigma) = 0;

                % build the derivative in x,y,z directions
                dK.x = obj.build_derivative(obj.lh_grid.V(:,1).', obj.rh_grid.V(:,1).', K, lh_dK, rh_dK, K_norm);
                dK.y = obj.build_derivative(obj.lh_grid.V(:,2).', obj.rh_grid.V(:,2).', K, lh_dK, rh_dK, K_norm);
                dK.z = obj.build_derivative(obj.lh_grid.V(:,3).', obj.rh_grid.V(:,3).', K, lh_dK, rh_dK, K_norm);
            end
        end
    end

    methods (Access = private)
        function dK_p = build_derivative(obj, lh_p, rh_p, K, lh_dK, rh_dK, K_norm)
            % get the directional derivative dK_p
            dK_p = blkdiag(lh_dK .* lh_p, rh_dK .* rh_p);

            % propagate the normalisation to the derivative
            dK_p = (dK_p - K .* sum(dK_p, 2)) ./ K_norm;

            % move to the GPU 
            dK_p = obj.maybe_GPU(dK_p);
        end

        function X = calc_distance(~, X)
            % calculate the distance (in mm) on an idealised brain.
            X = acos(max(min(X, 1), -1)) * 67;
        end
    end
end

