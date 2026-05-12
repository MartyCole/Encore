classdef QuarticKernel < Kernel
    % QuarticKernel
    %
    % Implements the quartic (biweight) kernel on the unit sphere S^2. This
    % avoids the singular derivative of acos(p·q) by mapping each neighbor 
    % q into the tangent plane at p using the log map: 
    %
    %       v = log_p(q) = theta * u
    %
    % where u is the unit tangent direction and theta is the geodesic angle.
    %
    % This class provides:
    %   - fast evaluation of the quartic kernel matrix K for a given bandwidth σ
    %   - directional derivatives of the kernel (x, y, z components)
    %
    % Internally, the kernel is computed using:
    %       K = (15/16σ) {1 - (r/σ)^2}^2 for r ≤ σ
    %
    % This class subclasses Kernel and implements compute_sigma.    

    properties (Access = private)
        lh_r
        lh_u
        rh_r
        rh_u
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

            % cache the distances vertices
            lh_dot = obj.clamp_corr(lh_grid.V * lh_grid.V.');
            rh_dot = obj.clamp_corr(rh_grid.V * rh_grid.V.');

            [obj.lh_r, obj.lh_u] = obj.logmap(lh_grid.V, lh_dot);
            [obj.rh_r, obj.rh_u] = obj.logmap(rh_grid.V, rh_dot);
        end

        function [K, dK] = compute_sigma(obj, sigma)
            % Compute the quartic kernel matrix K for bandwidth sigma.        
            %
            % Outputs:
            %   K  : kernel matrix
            %   dK : struct of directional derivatives (if requested)   
            
            % compute the kernel with the given sigma
            lh_K = (15/16) * (1/sigma) * (1 - (obj.lh_r ./ sigma).^2).^2;
            rh_K = (15/16) * (1/sigma) * (1 - (obj.rh_r ./ sigma).^2).^2;

            % only support within sigma
            lh_K(obj.lh_r > sigma) = 0;
            rh_K(obj.rh_r > sigma) = 0;

            K = blkdiag(lh_K, rh_K);

            % row normalise the kernel (diffusion smoothing)
            K_norm = sum(K,2);   
            K = K ./ K_norm;

            % move to the GPU?
            K = obj.maybe_GPU(K);

            % compute the directional derivative if requested
            if nargout == 2
                % compute the derivative of the kernel with the given sigma
                lh_dK = -(15/4) * (obj.lh_r ./ sigma.^3) .* (1 - (obj.lh_r ./ sigma).^2);
                rh_dK = -(15/4) * (obj.rh_r ./ sigma.^3) .* (1 - (obj.rh_r ./ sigma).^2);                       

                % only support within sigma
                lh_dK(obj.lh_r > sigma) = 0;
                rh_dK(obj.rh_r > sigma) = 0;

                % build the derivative in x,y,z directions
                dK.x = obj.build_derivative(lh_dK, rh_dK, obj.lh_u(:,:,1), obj.rh_u(:,:,1), K, K_norm);
                dK.y = obj.build_derivative(lh_dK, rh_dK, obj.lh_u(:,:,2), obj.rh_u(:,:,2), K, K_norm);
                dK.z = obj.build_derivative(lh_dK, rh_dK, obj.lh_u(:,:,3), obj.rh_u(:,:,3), K, K_norm);
            end
        end
    end

    methods (Access = private)
        function [r, u] = logmap(~, V, dotp)    
            % calculate angle (in mm) on an idealised brain.
            r = 67*acos(dotp); 
            
            V_expanded = reshape(V, [size(V,1), 1, 3]);           
            dotp_expanded = reshape(dotp, [size(dotp,1), size(dotp,2), 1]);

            % projection of q onto tangent plane at p
            proj = dotp_expanded .* V_expanded; 

            Vj = reshape(V, [1, size(V,1), 3]);  
            diff = Vj - proj;                   

            % normalize to get tangent direction
            normdiff = sqrt(sum(diff.^2, 3));   
            u = diff ./ max(normdiff, 1e-12);
        end       

        function dK_p = build_derivative(obj, lh_dK, rh_dK, lh_p, rh_p, K, K_norm)
            % get the directional derivative dK_p
            dK_p = blkdiag(-(lh_dK .* lh_p), -(rh_dK .* rh_p));

            % propagate the normalisation to the derivative
            dK_p = (dK_p - K .* sum(dK_p, 2)) ./ K_norm;

            % move to the GPU 
            dK_p = obj.maybe_GPU(dK_p);
        end

        function X = clamp_corr(~, X)
            X = max(min(X, 1), -1);
        end
    end
end
