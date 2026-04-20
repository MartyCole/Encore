classdef SphericalHeatKernel < Kernel
    % SphericalHeatKernel
    %
    % Implements the spherical heat kernel on the unit sphere S^2 using a
    % truncated spherical harmonic expansion up to degree H.
    %
    % This class provides:
    %   - fast evaluation of the heat kernel matrix K for a given bandwidth σ
    %   - directional derivatives of the kernel (x, y, z components)
    %   - sparse kernel density estimation (KDE) between point sets
    %   - efficient cross‑validation via cached Legendre polynomials
    %
    % Internally, the kernel is computed using:
    %       K(x,y) = Σ_{l=0}^H (2l+1) e^{-l(l+1)σ} P_l(x·y)
    %
    % Cached quantities:
    %   - Legendre polynomials P_l(cosθ) and derivatives dP_l/dx
    %   - cosine similarity matrices for LH and RH grids
    %   - precomputed Legendre tables for cutoff estimation
    %
    % This class subclasses Kernel and implements compute_sigma.

    properties (SetAccess = private)
        H
    end

    properties (Access = private)
        lh_X
        rh_X
        sigma_X

        % Legendre Polynomials
        lh_P
        lh_dP
        rh_P
        rh_dP

        sigma_P
        sigma_dP
    end

    methods
        function obj = SphericalHeatKernel(lh_grid, rh_grid, H, gpu)
            % Construct a spherical heat kernel object.
            %
            % Inputs:
            %   lh_grid, rh_grid : SphericalGrid objects for LH and RH
            %   H                : maximum harmonic degree for kernel expansion
            %   gpu              : (true) whether to store kernel matrices on GPU                   

            if nargin < 4
                gpu = true;
            end

            obj.use_GPU = gpu;

            obj.lh_grid = lh_grid;
            obj.rh_grid = rh_grid;
            obj.H = H;

            % cache the cross harmonic terms on the given grids for fast CV
            obj.lh_X = obj.clamp_corr(lh_grid.V * lh_grid.V.');
            obj.rh_X = obj.clamp_corr(rh_grid.V * rh_grid.V.');

            [obj.lh_P,obj.lh_dP] = obj.compute_legendre(obj.lh_X);
            [obj.rh_P,obj.rh_dP] = obj.compute_legendre(obj.rh_X);

            % cache Legendre polynomials for quick cutoff calculation
            obj.sigma_X = linspace(1, -1, 5000);
            [obj.sigma_P, ~] = obj.compute_legendre(obj.sigma_X);
        end      

        function [K, dK] = compute_sigma(obj, sigma)
            % Compute the heat kernel matrix K for bandwidth sigma.        
            %
            % Outputs:
            %   K  : kernel matrix
            %   dK : struct of directional derivatives (if requested)

            cutoff_dst = obj.compute_cutoff(sigma);

            % compute the kernel
            lh_K = obj.calc_cross_harmonics(obj.lh_P, sigma, cutoff_dst, obj.lh_X);
            rh_K = obj.calc_cross_harmonics(obj.rh_P, sigma, cutoff_dst, obj.rh_X);

            K = blkdiag(lh_K, rh_K);

            % row normalise the kernel (diffusion smoothing)
            K_norm = sum(K,2);   
            K = K ./ K_norm;

            % move to the GPU?
            K = obj.maybe_GPU(K);

            % compute the directional derivative if requested
            if nargout == 2
                % note that here we use the derivatives of Legendre polynomials
                lh_dKp = obj.calc_cross_harmonics(obj.lh_dP, sigma, cutoff_dst, obj.lh_X);
                rh_dKp = obj.calc_cross_harmonics(obj.rh_dP, sigma, cutoff_dst, obj.rh_X);

                % build the derivative in x,y,z directions
                dK.x = obj.build_derivative(lh_dKp.', rh_dKp, obj.lh_grid.V(:,1).', obj.rh_grid.V(:,1).', K, K_norm);
                dK.y = obj.build_derivative(lh_dKp.', rh_dKp, obj.lh_grid.V(:,2).', obj.rh_grid.V(:,2).', K, K_norm);
                dK.z = obj.build_derivative(lh_dKp.', rh_dKp, obj.lh_grid.V(:,3).', obj.rh_grid.V(:,3).', K, K_norm);
            end
        end
    end

    methods (Access = private)
        function cutoff_dst = compute_cutoff(obj, sigma)
            % Estimate cosine cutoff threshold for kernel truncation.
            %
            % Returns:
            %   cutoff_dst : smallest cosine value for which kernel > epsilon

            Hs = 0:obj.H;
            epsilon = 1e-3;

            % harmonic weights
            L = Hs .* (Hs + 1);
            W0 = ((2*Hs+1)' * (2*Hs+1));

            W = exp(-sigma * (L + L')) .* W0;
            kernel_vals = obj.sigma_P.' * W * obj.sigma_P(:,1);

            idx = find(kernel_vals < epsilon, 1, 'first');            
            idx = find(kernel_vals(1:idx) > 0, 1, 'last');

            if isempty(idx) 
                cutoff_dst = -1;
            else
                cutoff_dst = obj.sigma_X(idx);
            end
        end

        function [K] = calc_cross_harmonics(obj, P, sigma, cutoff_dst, X)
            % Evaluate truncated spherical harmonic expansion for all pairs.
            %
            % Inputs:
            %   P          : Legendre polynomials P_l(x·y)
            %   sigma      : bandwidth
            %   cutoff_dst : cosine cutoff threshold
            %   X          : cosine similarity matrix
            %
            % Output:
            %   K : kernel values with cutoff applied

            Hs = [0:obj.H].';

            K = squeeze(sum(((2*Hs)+1) .* exp(-Hs.*(Hs+1)*sigma) .* P, 1));

            K(X < cutoff_dst) = 0;
        end

        function dK_p = build_derivative(obj, lh_dK, rh_dK, lh_p, rh_p, K, K_norm)
            % Construct directional derivative of the kernel in Cartesian directions.
            %
            % Inputs:
            %   lh_dK, rh_dK : derivative of Legendre expansion for LH and RH
            %   lh_p, rh_p   : x, y, or z coordinates of grid vertices
            %   K            : the Kernel to take the derivative of
            %   K_norm       : the constant used to normalise K
            %
            % Output:
            %   dK_p : block‑diagonal derivative matrix

            % get the directional derivative dK_p
            dK_p = blkdiag(lh_dK .* lh_p, rh_dK .* rh_p);

            % propagate the normalisation to the derivative
            dK_p = (dK_p - K .* sum(dK_p, 2)) ./ K_norm;

            % move to the GPU 
            dK_p = obj.maybe_GPU(dK_p);
        end
    
        function [P, dP] = compute_legendre(obj, X)
            % Compute Legendre polynomials P_l(x) and derivatives dP_l/dx
            % for all l = 0..H using stable recurrence relations.
            %
            % Input:
            %   X : cosine values (scalar, vector, or matrix)
            %
            % Outputs:
            %   P  : (H+1)×size(X) array of Legendre values
            %   dP : (H+1)×size(X) array of derivatives

            N = size(X);

            if N(1) == 1
                N(1) = N(2);
                N(2) = 1;
            end

            P = zeros([obj.H+1, N]);
            dP = zeros([obj.H+1, N]);

            P(1,:,:) = 1;
            dP(1,:,:) = 0;

            if obj.H >= 1
                P(2,:,:) = X;
                dP(2,:,:) = 1;
            end

            % initialise the recurrence
            P_hm1 = squeeze(P(1,:,:));
            P_h = squeeze(P(2,:,:));
            dP_hm1 = squeeze(dP(1,:,:));
            dP_h = squeeze(dP(2,:,:));

            for h = 1:obj.H-1
                % recurrence for P_{h+1} 
                P_hp1 = ((2*h+1).*X.*P_h - h.*P_hm1) ./ (h+1);
            
                % derivative recurrence (version without division of 1-x^2) 
                dP_hp1 = ((2*h+1).*(P_h + X.*dP_h) - h.*dP_hm1) ./ (h+1);
            
                % store results
                P(h+2,:,:) = P_hp1;
                dP(h+2,:,:) = dP_hp1;
            
                % shift recurrence for next step
                P_hm1 = P_h;
                P_h = P_hp1;
                dP_hm1 = dP_h;
                dP_h = dP_hp1;
            end
        end     

        function X = clamp_corr(~, X)
            % Clamp correlation (cosine similarity) values to the range [-1, 1].
            %
            % Prevents numerical issues in Legendre evaluation.
            X = max(min(X, 1), -1);
        end
    end
end

