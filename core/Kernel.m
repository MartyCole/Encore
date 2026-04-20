classdef (Abstract) Kernel < handle
    % Kernel (Abstract)
    %
    % This abstract class defines the interface and shared functionality for
    % spherical smoothing kernels used in connectome estimation and registration.
    %
    % A Kernel object is responsible for:
    %   - computing a kernel matrix K for a given bandwidth sigma
    %   - optionally computing kernel derivatives dK (x, y, z directions)
    %   - supporting GPU acceleration when enabled
    %   - performing cross‑validation to estimate the optimal bandwidth
    %
    % Subclasses must implement:
    %       [K, dK] = compute_sigma(obj, sigma)
    %
    % Key Methods:
    %   cross_validate_sigma : estimate optimal sigma via ISE minimization
    %   compute_ISE          : compute leave‑one‑out ISE for a kernel
    %
    % This class provides the common infrastructure for all kernel types
    % (e.g., heat kernels, von Mises–Fisher kernels, diffusion kernels),
    % while allowing subclasses to define the specific kernel formulation.
    %
    properties
        lh_grid % spherical grid for the left hemisphere
        rh_grid % spherical grid for the right hemisphere

        use_GPU logical = true % convert matrices to GPU when enabled
    end

    methods (Abstract)
        % Abstract method to compute the kernel matrix and its derivatives.
        %
        % Inputs:
        %   sigma : kernel bandwidth parameter
        %
        % Outputs:
        %   K  : kernel matrix evaluated at bandwidth sigma
        %   dK : struct containing kernel derivatives (x, y, z)
        %
        [K, dK] = compute_sigma(obj, sigma)        
    end

    methods
        function [sigma_opt, LL] = cross_validate_sigma(obj, concon, sigmas)
            % Perform cross‑validation to estimate the optimal kernel bandwidth sigma.
            %
            % Inputs:
            %   concon : Concon object containing fiber endpoints
            %   sigmas : list of candidate sigma values
            %
            % Outputs:
            %   sigma_opt : sigma value minimizing the integrated squared error (ISE)
            %   LL        : Leave one out Log Liklihood values for each candidate sigma       

            LL = zeros(size(sigmas));

            % because lots of small multiplications, quicker in CPU
            A = gather(concon.get_adjacency());

            fprintf("-------------------------------------------\n");
            fprintf("Performing bandwidth estimation:\n");           
            fprintf("-------------------------------------------\n");

            S = numel(sigmas);

            % for printing the progress
            thr = (S > 20) * round(S / 20) + (S <= 20); 
       
            for i = 1:numel(sigmas)
                sigma = sigmas(i);

                % because lots of small multiplications, quicker in CPU
                K = gather(obj.compute_sigma(sigma));

                % calculate the LL for the current sigma
                LL(i) = obj.compute_LOO(K, A, concon);

                % print progress
                if mod(i/S,thr) == 0
                    fprintf("=")
                end
            end

            % return the sigma with maximum logliklihood
            [~, idx] = max(LL);
            sigma_opt = sigmas(idx);
            
            fprintf("\nOptimal sigma = %f\n", sigma_opt)           
        end
    end

    methods (Access = protected)
        function LOO_val = compute_LOO(~, K, A, concon)
            % Compute the Leave-one-out Log-Liklihood for a given kernel K.
            %
            % Inputs:
            %   K          : kernel matrix
            %   A          : adjacency matrix
            %   concon     : Concon object with barycentric endpoint data
            %   chunk_size : optional batch size for memory‑efficient processing
            %
            % Output:
            %   LOO_val    : mean LOO-LL value across fibers
            
            % number of fibers
            M = size(concon.st_points,1);

            % full estimated connectome
            C = max(K.' * A * K, 0);
          
            % save local for speed
            W_sp = concon.st_points;
            I_sp = concon.st_idx;
            W_ep = concon.en_points;
            I_ep = concon.en_idx;

            % somewhere to put the results
            vals = zeros(M,1);

            for s = 1:M
                % fiber coordinates and indices
                t_sp = I_sp(s,:); t_ep = I_ep(s,:);   
                w_sp = W_sp(s,:); w_ep = W_ep(s,:);   
              
                % the fiber's contribution to the kernel
                K_local = (K(t_sp, t_sp) * K(t_ep, t_ep).') / (M - 1);

                % subtract the fiber from the kernel and interpolate
                vals(s) = max(w_sp * (C(t_sp,t_ep) - K_local) * w_ep.', eps);
            end
 
            % final leave-one-out log liklihood
            LOO_val = mean(log(vals));
        end

        function K = maybe_GPU(obj, K)
            % Convert kernel matrix to GPU if enabled.
            %
            % Ensures:
            %   - single precision
            %   - gpuArray conversion only when use_GPU is true
            %
            % Output:
            %   K : possibly GPU‑resident kernel matrix

            K = single(K);

            if obj.use_GPU == true
                K = gpuArray(K);
            end
        end 
    end
end

