classdef GaussianKernel < Kernel
    % GaussianKernel
    %
    % Implements the Gaussian kernel on the unit sphere S^2
    %
    % This class provides:
    %   - fast evaluation of the quartic kernel matrix K for a given bandwidth σ
    %   - directional derivatives of the kernel (x, y, z components)
    %
    % Internally, the kernel is computed using:
    %       K = (1/(2πσ^2)) e^{-(X^2 / (2σ^2})}
    %
    % This class subclasses Kernel and implements compute_sigma.

    properties (Access = private)
        lh_X
        rh_X    

        epsilon
    end

    methods
        function obj = GaussianKernel(lh_grid, rh_grid, epsilon, gpu)
            % Construct a Gaussian kernel object.
            %
            % Inputs:
            %   lh_grid, rh_grid : SphericalGrid objects for LH and RH          
            %   gpu              : (true) whether to store kernel matrices on GPU                   

            if nargin < 4
                gpu = true;
            end

            obj.use_GPU = gpu;
            obj.epsilon = epsilon;

            % save the grid the kernel is on
            obj.lh_grid = lh_grid;
            obj.rh_grid = rh_grid;    

            % cache the distance between vertices on the grid
            obj.lh_X = obj.calc_distance(lh_grid.V * lh_grid.V.');
            obj.rh_X = obj.calc_distance(rh_grid.V * rh_grid.V.');
        end      

        function [K, dK] = compute_sigma(obj, sigma)
            % Compute the Gaussian kernel matrix K for bandwidth sigma.        
            %
            % Outputs:
            %   K  : kernel matrix
            %   dK : struct of directional derivatives (if requested)    

            % calculate cutoff distance based on FWHM and epsilon
            cutoff_dst = (sigma * sqrt(8*log(2))) * sqrt((-log2(obj.epsilon)));
           
            % compute the kernel with the given sigma
            lh_K = (1/(2*pi*sigma^2)) .* exp(-(obj.lh_X.^2 / (2 * (sigma^2))));
            rh_K = (1/(2*pi*sigma^2)) .* exp(-(obj.rh_X.^2 / (2 * (sigma^2))));

            % only support within cutoff distance
            lh_K(obj.lh_X > cutoff_dst) = 0;
            rh_K(obj.rh_X > cutoff_dst) = 0;
            
            K = blkdiag(lh_K, rh_K);            

            % row normalise the kernel (diffusion smoothing)
            K_norm = sum(K,2);   
            K = K ./ K_norm;

            % move to the GPU?
            K = obj.maybe_GPU(K);

            % compute the directional derivative if requested
            if nargout == 2      
                % compute the derivative of the kernel with the given sigma
                lh_dK = (-obj.lh_X/(2*pi*sigma^4)) .* exp(-(obj.lh_X.^2 / (2 * (sigma^2))));
                rh_dK = (-obj.rh_X/(2*pi*sigma^4)) .* exp(-(obj.rh_X.^2 / (2 * (sigma^2))));

                % only support within cutoff distance
                lh_dK(obj.lh_X > cutoff_dst) = 0;
                rh_dK(obj.rh_X > cutoff_dst) = 0;

                % build the derivative in x,y,z directions
                dK.x = obj.build_derivative(obj.lh_grid.V(:,1).', obj.rh_grid.V(:,1).', K, lh_dK, rh_dK, K_norm);
                dK.y = obj.build_derivative(obj.lh_grid.V(:,2).', obj.rh_grid.V(:,2).', K, lh_dK, rh_dK, K_norm);
                dK.z = obj.build_derivative(obj.lh_grid.V(:,3).', obj.rh_grid.V(:,3).', K, lh_dK, rh_dK, K_norm);
            end
        end

        function sigma = FWHM_to_sigma(~, FWHM)
            sigma = (FWHM/2) / sqrt(8*log(2));
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

