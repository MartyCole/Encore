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
        function [sigma_opt, ISE] = cross_validate_sigma(obj, concon, sigmas)
            % Perform cross‑validation to estimate the optimal kernel bandwidth sigma.
            %
            % Inputs:
            %   concon : Concon object containing fiber endpoints
            %   sigmas : list of candidate sigma values
            %
            % Outputs:
            %   sigma_opt : sigma value minimizing the integrated squared error (ISE)
            %   ISE       : ISE values for each candidate sigma       

            ISE = zeros(size(sigmas));

            A = gather(concon.get_adjacency());

            fprintf("-------------------------------------------\n");
            fprintf("Performing bandwidth estimation:\n");           
            fprintf("-------------------------------------------\n");

            for i = 1:numel(sigmas)
                sigma = sigmas(i);

                K = gather(obj.compute_sigma(sigma));

                ISE(i) = obj.compute_ISE(K, A, concon);

                if mod(i,(numel(sigmas)/20)) == 0
                    fprintf("=")
                end
            end

            [~, idx] = min(ISE);
            sigma_opt = sigmas(idx);
            
            fprintf("\nOptimal sigma = %f\n", sigma_opt)           
        end
    end

    methods (Access = protected)
        function ISE_val = compute_ISE(~, K, A, concon, chunk_size)
            % Compute the Integrated Squared Error (ISE) for a given kernel K.
            %
            % Inputs:
            %   K          : kernel matrix
            %   A          : adjacency matrix
            %   concon     : Concon object with barycentric endpoint data
            %   chunk_size : optional batch size for memory‑efficient processing
            %
            % Output:
            %   ISE_val : mean ISE across all fibers
            
            if nargin < 5
                chunk_size = 5000;
            end

            M = size(concon.st_points,1);

            C = max(K.' * A * K, 0);
           
            vals = zeros(M,1);

            W_sp = concon.st_points;
            I_sp = concon.st_idx;
            W_ep = concon.en_points;
            I_ep = concon.en_idx;

            % process pairs in chunks
            for s0 = 1:chunk_size:M
                s1  = min(s0 + chunk_size - 1, M);
                idx = s0:s1;
                B   = numel(idx);
             
                for b = 1:B
                    s = idx(b);

                    t_sp = I_sp(s,:); 
                    t_ep = I_ep(s,:);   

                    w_sp = W_sp(s,:);   
                    w_ep = W_ep(s,:);   

                    full_block = C(t_sp, t_ep);                    
                    self_block = (K(t_sp,t_sp)' * ((w_sp'*w_ep)/(2*M)) * K(t_ep,t_ep));                                                  
           
                    vals(s) = (w_sp * full_block.^2 * w_ep.').^2 - 2*(w_sp * (full_block - self_block) * w_ep.');  
                end
            end

            ISE_val = mean(vals);
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

