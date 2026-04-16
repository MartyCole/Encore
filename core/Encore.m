classdef Encore
    % Encore    mEtric-based coNtinuous-COnnectivity REgistration
    %
    % This class performs spherical registration between two structural
    % connectomes (or general spherical functions). It supports:
    %
    %   - computing a Karcher mean template across subjects
    %   - performing linear registration on the sphere (rotation)
    %   - performing non‑linear diffeomorphic registration on the sphere
    %
    % The class operates on two spherical grids (LH and RH), computes
    % integration weights, and iteratively aligns a moving function F2 to a
    % fixed function F1 using a kernel representation.
    %
    properties (GetAccess = private)
        lh_grid % Left hemisphere spherical grid
        rh_grid % Right hemisphere spherical grid
        A % PxP integration matrix (triangle areas)
    end

    properties (Access = public)        
        delta % Step size for gradient descent
        max_iters % Maximum number of registration iterations
        threshold % Convergence threshold on cost reduction    

        use_GPU logical = true % convert matrices to GPU when enabled
    end

    methods
        function obj = Encore(lh_mesh,rh_mesh,l,delta,iters,threshold,gpu)
            % Constructor for Encore
            %
            % Inputs:
            %   lh_mesh, rh_mesh : spherical meshes for LH and RH
            %   l                : spherical harmonic bandwidth / resolution
            %   delta            : gradient descent step size
            %   iters            : maximum iterations
            %   threshold        : convergence threshold
            %
            % Builds spherical grids and precomputes the integration matrix A.

            if nargin < 7
                gpu = true;
            end
            
            obj.use_GPU = gpu;

            obj.lh_grid = SphericalGrid(lh_mesh,l,gpu);
            obj.rh_grid = SphericalGrid(rh_mesh,l,gpu);

            obj.A = obj.maybe_GPU([obj.lh_grid.A; obj.rh_grid.A] * [obj.lh_grid.A; obj.rh_grid.A].');

            obj.delta = delta;
            obj.max_iters = iters;
            obj.threshold = threshold;            
        end

        function template = get_template(~, Fs, kernel, iters)
            % Compute a Karcher mean template from a set of subjects.
            %
            % Inputs:
            %   Fs     : A cell array containing a population of Concon objects
            %   kernel : The kernel to evaluated connectomes with
            %   iters  : The maximum number of iterations of the algorithm            
            %
            % Output:
            %   template : the estimated mean Q function

            N_subs = size(Fs,1);
                        
            % find the mean Q function            
            for i = 1:N_subs
                tmp = sqrt(Fs{i}.evaluate(kernel));                

                if i == 1
                    Q_bar = tmp;
                else
                    Q_bar = Q_bar + tmp;
                end
                disp(i)
            end
            
            Q_bar = Q_bar ./ N_subs;
            Q_norm = zeros(1,N_subs);
            
            % find the Q function nearest to the mean
            for i = 1:N_subs
                Q_norm(i) = sum((sqrt(Fs{i}.evaluate(kernel)) - Q_bar).^2, 'all');
                disp(i)
            end
            
            idx = find(Q_norm == min(Q_norm),1);
            Q_mu = sqrt(Fs{idx}.evaluate(kernel));
            
            % iterate closer the the karcher median
            for iter = 1:iters
                vv = zeros(size(Q_mu,1),size(Q_mu,1),size(Fs,1));
                distance = zeros(N_subs,1);
                    
                for i = 1:N_subs
                    tmpQ = sqrt(Fs{i}.evaluate(kernel));
                    
                    tmpTheta = sum(tmpQ(:) .* Q_mu(:));
                    
                    if 1 - abs(tmpTheta) < 1e-14
                        tmpTheta = sign(tmpTheta);
                    end

                    distance(i) = acos(tmpTheta);
                    
                    if distance(i) > 0
                        vv(:,:,i) = ((distance(i) / sin(distance(i))) * (tmpQ - cos(distance(i))*Q_mu)) / distance(i);
                    end
                end
                
                v_bar = sum(vv,3) / sum(1./distance(distance > 0));
                tmp = sqrt(sum(v_bar(:) .* v_bar(:)));
                Q_mu = (cos(0.2*tmp) * Q_mu) + (sin(0.2*tmp) * (v_bar / tmp));
                Q_mu = Q_mu / sqrt(sum(Q_mu(:) .* Q_mu(:)));

                fprintf('Template norm: %0.4f\n', tmp);

                if tmp < 0.005
                    break
                end
            end
            
            template = Q_mu;
        end

        function [lh_warp,rh_warp,cost] = register(obj,F1,F2,kernel,kernel_diff,varargin)         
            % Register F2 (moving) to F1 (fixed) using non‑linear spherical warping.
            %
            % Inputs:
            %   F1, F2       : Concon objects or precomputed Q functions
            %   kernel       : smoothing kernel
            %   kernel_diff  : kernel derivatives (x,y,z)
            %
            % Optional:
            %   'verbose'        : if greater than 0, plot progress
            %   'init_rotation'  : whether to perform initial rigid alignment
            %
            % Outputs:
            %   lh_warp, rh_warp : resulting spherical warp fields
            %   cost             : cost function values per iteration            

            p = inputParser;
            addParameter(p, 'verbose', 0, @isnumeric);
            addParameter(p, 'init_rotation', false, @islogical);

            % parse optional variables
            parse(p, varargin{:});
            params = p.Results;

            % setup warps and variables for the loop
            F2 = copy(F2);

            if params.init_rotation == true
                [lh_warp, rh_warp] = obj.initial_registration(F1,F2,kernel);

                F2.warp_connectome(lh_warp, rh_warp); 
            else
                lh_warp = SphericalWarp(obj.lh_grid);
                rh_warp = SphericalWarp(obj.rh_grid);    
            end            

            lh_P = length(obj.lh_grid.V);
            all_P = length(obj.lh_grid.V) + length(obj.rh_grid.V);

            idx_a = 1:lh_P; 
            idx_b = (lh_P+1):all_P;

            lh_idx = [idx_a,idx_b];
            rh_idx = [idx_b,idx_a];

            % initial FIXED function     
            if isa(F1,"Concon")
                Q1 = sqrt(F1.evaluate(kernel));                
            else                
                Q1 = F1;
            end
            
            % initial MOVING function
            [Q2,Q2_e1,Q2_e2] = F2.evaluate(kernel, kernel_diff);
            [Q2,Q2_e1,Q2_e2] = obj.Q_transform(Q2,Q2_e1,Q2_e2);     
           
            % initial cost
            FmM = (Q1 - Q2);
            cost = zeros(obj.max_iters+1,1);   

            cost(1) = sum(FmM.^2,'all');

            fprintf("-------------------------------------------\n");
            fprintf("Performing non-linear registration:\n");
            fprintf('Initial Cost: %0.6f\n', cost(1))           
            
            % start registering
            for iter = 1:obj.max_iters  
                % compute left and right gradients
                lh_dH = obj.compute_gradient(FmM, Q2, Q2_e1, Q2_e2, idx_a, lh_idx, obj.lh_grid);               
                rh_dH = obj.compute_gradient(FmM, Q2, Q2_e1, Q2_e2, idx_b, rh_idx, obj.rh_grid);

                % constrain step size to avoid folded triangles
                lh_step_size = obj.compute_step_size(lh_dH, obj.delta, 0.2);
                rh_step_size = obj.compute_step_size(rh_dH, obj.delta, 0.2);                

                % update the warps
                lh_warp.compose_warp(lh_step_size);                 
                rh_warp.compose_warp(rh_step_size);  

                % re-evaluate Q2 and its derivatives
                F2.warp_connectome(lh_warp, rh_warp);                
               
                [Q2,Q2_e1,Q2_e2] = F2.evaluate(kernel, kernel_diff);
                [Q2,Q2_e1,Q2_e2] = obj.Q_transform(Q2,Q2_e1,Q2_e2);
                
                % cost update
                FmM = (Q1 - Q2); 
                cost(iter+1) = sum(FmM.^2,'all');             
                
                % check ending condition
                if ((cost(iter) - cost(iter+1)) < obj.threshold)       
                   fprintf('Converged after %d iterations\n', iter)            
                   break
                end                                  

                if (params.verbose > 0)
                    obj.plot_progress(params.verbose, lh_warp, rh_warp, cost, iter);                  
                elseif mod(iter,10) == 0
                    fprintf("Current Cost: %0.6f - %i iterations\n", cost(iter), iter)
                end
            end                      

            fprintf("Final Cost:   %0.6f\n", cost(iter)); 
            fprintf("-------------------------------------------\n");
           
            % evaluate function after warping
            cost = cost(1:iter+1);
        end
    end

    methods (Access = protected)
        function [lh_warp, rh_warp] = initial_registration(obj,F1,F2,kernel)
            % Perform rigid (rotation‑only) registration between two connectomes
            % before non‑linear warping begins.
            %
            % This step finds the best pair of SO(3) rotations for LH and RH that
            % minimize the squared distance between Q1 and Q2. The search uses:
            %   - icosahedral permutations (60 rotational symmetries)
            %   - multi‑shell incremental rotations (precomputed)
            %   - barycentric interpolation for rotated matrices
            %
            % Inputs:
            %   F1     : The fixed Concon object
            %   F2     : The moving Concon object
            %   kernel : The kernel to evaluated connectomes with
            %
            % Output:
            %   lh_warp, rh_warp : SphericalWarp objects containing the best rigid
            %                      rotation for each hemisphere.

            % initial FIXED function     
            if isa(F1,"Concon")
                Q1 = sqrt(F1.evaluate(kernel));                
            else                
                Q1 = F1;
            end

            %Q1 = sqrt(F1.evaluate(kernel));
            Q2 = sqrt(F2.evaluate(kernel));

            fprintf("-------------------------------------------\n");
            fprintf("Performing rigid registration:\n");
            fprintf("Initial Cost: %0.6f\n", sum((Q1-Q2).^2, 'all'));   
            fprintf("-------------------------------------------\n");

            % load orthogonal rotations in SO(3)
            lh_perm = icosahedral_permutations(obj.lh_grid.V);
            rh_perm = icosahedral_permutations(obj.rh_grid.V);
            ico_rot = icosahedral_rotations();
          
            % save important indices
            lh_idx = 1:size(obj.lh_grid.V,1);
            rh_idx = size(obj.lh_grid.V,1) + (1:size(obj.rh_grid.V,1));

            % we deal with left and right hemispheres individually
            lh_Q1 = Q1(lh_idx, lh_idx);
            rh_Q1 = Q1(rh_idx, rh_idx);

            % load precomputed shells
            classFolder = fileparts(mfilename('fullpath'));
            load(fullfile(classFolder, 'data/rigid_search_grid.mat'), 'shells');

            % ---------------------------------------------------
            % Search shell 1 combined with the 60 SO(3) rotations
            % ---------------------------------------------------
            fprintf("Searching Shell: %i -> ", 1);

            R = shells{1}.R;   
            
            lh_costs = zeros(60, numel(R));
            rh_costs = zeros(60, numel(R));
            
            for s = 1:numel(R)
                % interpolate rotations in the shell
                Q2_rot = F2.rotate_matrix_barycentric(Q2, R{s}, R{s});

                % permute through the SO(3) rotations
                lh_Q2 = Q2_rot(lh_idx, lh_idx);
                rh_Q2 = Q2_rot(rh_idx, rh_idx);
                            
                for i = 1:60                    
                    lh_costs(i,s) = sum((lh_Q1 - lh_Q2(lh_perm(:,i),lh_perm(:,i))).^2,'all');
                    rh_costs(i,s) = sum((rh_Q1 - rh_Q2(rh_perm(:,i),rh_perm(:,i))).^2,'all');
                end

                % print progress
                if mod(s,6) == 0
                    fprintf("=");
                end
            end
            
            fprintf("\n");
           
            % select the best shells{1}.K candidate rotations to continue the search from
            [lh_candidate, rh_candidate] = obj.select_best_candidates(lh_costs, rh_costs, shells{1}.K, R, ico_rot);           
            
            % ------------------------------------------------
            % Search remaining shells using just interpolation
            % ------------------------------------------------
            for sh = 2:numel(shells)
                fprintf("Searching Shell: %i -> ", sh);

                R = shells{sh}.R;
                prev_K = shells{sh-1}.K;
                curr_K = shells{sh}.K;
            
                lh_costs = zeros(prev_K, numel(R));
                rh_costs = zeros(prev_K, numel(R));
            
                for k = 1:prev_K
                    for s = 1:numel(R)
                        % interpolate rotations in the shell composed with the candidate intital rotation
                        Q2_rot = F2.rotate_matrix_barycentric(Q2, lh_candidate{k} * R{s}, rh_candidate{k} * R{s});
            
                        lh_costs(k,s) = sum((lh_Q1 - Q2_rot(lh_idx, lh_idx)).^2, 'all');
                        rh_costs(k,s) = sum((rh_Q1 - Q2_rot(rh_idx, rh_idx)).^2, 'all');
                     
                        % print progress
                        if mod(s,prev_K) == 0
                            fprintf("=");
                        end                  
                    end
                end

                fprintf("\n");
            
                % select the best shells{1}.K candidate rotations to continue the search from
                [lh_candidate, rh_candidate] = obj.select_best_candidates(lh_costs, rh_costs, curr_K, R, lh_candidate, rh_candidate); 
            end

            lh_warp = SphericalWarp(obj.lh_grid);
            rh_warp = SphericalWarp(obj.rh_grid);   

            lh_warp.rotate(lh_candidate{1}');
            rh_warp.rotate(rh_candidate{1}');            
        end

        function [lh_new, rh_new] = select_best_candidates(~, lh_costs, rh_costs, K, R_list, varargin)
            % Select the top-K rotation candidates for the next shell of the
            % multi-shell rotation search.
            %
            % Supports two modes:
            %
            %   Shell 1:
            %       - Uses icosahedral rotations (ico_rot)
            %       - No previous candidates exist
            %
            %   Shells 2–5:
            %       - Refines previous candidates by composing them with
            %         rotations from R_list
            %
            % Inputs:
            %   lh_costs, rh_costs : cost matrices for LH and RH rotations
            %   K                  : number of best candidates to keep
            %   R_list             : list of incremental rotations
            %
            % Optional:
            %   Shell 1:
            %       varargin{1} = ico_rot (3×3×N rotation matrices)
            %
            %   Shells 2–5:
            %       varargin{1} = previous LH candidates
            %       varargin{2} = previous RH candidates
            %
            % Outputs:
            %   lh_new, rh_new : cell arrays of K new rotation candidates

            % determine if this is Shell 1 or later shells
            if nargin == 6
                % shell 1: use ico_rot
                ico_rot = varargin{1};
                lh_prev = [];
                rh_prev = [];
            else
                % shells 2–5: use previous candidates

                lh_prev = varargin{1};
                rh_prev = varargin{2};
            end

            % sort costs
            [~, lh_idx] = sort(lh_costs(:));
            [~, rh_idx] = sort(rh_costs(:));

            [lh_a, lh_b] = ind2sub(size(lh_costs), lh_idx(1:K));
            [rh_a, rh_b] = ind2sub(size(rh_costs), rh_idx(1:K));

            lh_new = cell(K,1);
            rh_new = cell(K,1);

            for k = 1:K
                if isempty(lh_prev)
                    % shell 1: combine R_list with icosahedral rotations
                    lh_new{k} = R_list{lh_b(k)} * ico_rot(:,:,lh_a(k))';
                    rh_new{k} = R_list{rh_b(k)} * ico_rot(:,:,rh_a(k))';
                else
                    % shells 2–5: refine previous candidates
                    lh_new{k} = lh_prev{lh_a(k)} * R_list{lh_b(k)};
                    rh_new{k} = rh_prev{rh_a(k)} * R_list{rh_b(k)};
                end
            end
        end

        function dH = compute_gradient(~, FmM, M, M_e1, M_e2, idx, hemi_idx, grid) 
            % Compute gradient of the registration cost with respect to the warp.
            %
            % Inputs:
            %   FmM      : (Q1 - Q2) difference
            %   M        : Q2
            %   M_e1,e2  : directional derivatives of Q2
            %   idx      : vertex indices for this hemisphere
            %   hemi_idx : indices of the opposite hemisphere block
            %   grid     : spherical grid (LH or RH)
            %
            % Output:
            %   dH       : gradient vector field on the sphere

            % extract the relevant hemisphere block
            FmM_loc  = FmM(idx, hemi_idx);
            
            % precompute the scalar contractions for each vertex            
            S1 = sum(FmM_loc .* M_e1(idx, hemi_idx), 2); 
            S2 = sum(FmM_loc .* M_e2(idx, hemi_idx), 2); 
            S3 = sum(FmM_loc .* M(idx, hemi_idx),    2);
            
            % build the integrand for each vertex
            integrand = 2 * S1 .* grid.basis(:, :, 1) ...
                      + 2 * S2 .* grid.basis(:, :, 2) + S3 .* grid.divergence;
                       
            % integrate over the sphere
            dH_vec = 2 * sum(integrand, 1);
            
            % project back onto the basis           
            dH = squeeze(sum(dH_vec .* grid.basis, 2));
        end

        function step_size = compute_step_size(~, dH, delta, threshold)
            % Compute a stable step size for gradient descent.
            %
            % Ensures the update does not fold triangles by limiting the maximum
            % vector magnitude.
            %
            % Output:
            %   step_size : scaled update field

            step_norm = max(vecnorm(dH,2,2));

            if step_norm > threshold
                dH = dH * (threshold / step_norm);
            end
        
            step_size = -delta*dH;           
        end

        function [Q, Q_e1, Q_e2] = Q_transform(~, F, F_e1, F_e2)
            % Convert F to Q = sqrt(F) and propagate derivatives via chain rule.

            Q = sqrt(F);

            % propegate the square root to the derivative by chain rule
            Q_e1 = ((1./(2*(max(Q,1e-15)))) .* F_e1);
            Q_e2 = ((1./(2*(max(Q,1e-15)))) .* F_e2);
        end

        function plot_progress(~, fig_id, lh_warp, rh_warp, cost, iter)
            % Visualize registration progress during optimization.
            %
            % Displays:
            %   - Estimated left-hemisphere warp field
            %   - Estimated right-hemisphere warp field
            %   - Cost function trajectory up to current iteration
            %
            % Inputs:
            %   fig_id   : figure handle or ID to reuse
            %   lh_warp  : SphericalWarp object for LH
            %   rh_warp  : SphericalWarp object for RH
            %   cost     : cost values across iterations
            %   iter     : current iteration index

            fig = figure(fig_id);
            t = tiledlayout(2,2);
            t.Padding = 'compact';
            t.TileSpacing = 'compact';

            nexttile
            lh_warp.plot('Estimated LH-Warp')
            clim([0 2])
            view([0 90 0]);
            nexttile
            rh_warp.plot('Estimated RH-Warp')
            clim([0 2])
            view([0 90 0]);

            nexttile([1,2])
            plot(1:(iter+1),cost(1:(iter+1)))
            xlim([1,100]);
            ylim([0,cost(1)]);
            title(sprintf('Iteration %d cost: %0.6f\n', iter, cost(iter+1)));

            fig.Units = "normalized";
            fig.OuterPosition = [0.154166666666667 0.08 0.411805555555556 0.834444444444444];
        end

        function X = maybe_GPU(obj, X)
            % Convert matrix to GPU if enabled.            
            %
            % Output:
            %   X : possibly GPU‑resident matrix

            if obj.use_GPU == true
                X = gpuArray(X);
            end
        end
    end
end
