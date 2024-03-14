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
classdef Encore
    properties (GetAccess = private)
        lh_grid
        rh_grid
        interpolator
        A
    end

    properties (Access = public)        
        delta
        max_iters
        threshold        
    end

    methods
        function obj = Encore(lh_mesh,rh_mesh,l,delta,iters,threshold)
            obj.lh_grid = SphericalGrid(lh_mesh,l);
            obj.rh_grid = SphericalGrid(rh_mesh,l);
             
            obj.interpolator = ConConInterpolator(obj.lh_grid, obj.rh_grid, 1e-10);

            obj.A = [obj.lh_grid.A; obj.rh_grid.A] * [obj.lh_grid.A; obj.rh_grid.A].';

            obj.delta = delta;
            obj.max_iters = iters;
            obj.threshold = threshold;            
        end

	function template = get_template(obj, Fs, iters)
            N_subs = size(Fs,3);
            
            % find the mean Q function            
            for i = 1:N_subs
                Fs(:,:,i) = sqrt(Fs(:,:,i) ./ sum(Fs(:,:,i) .* obj.A, 'all'));  
                disp(i)
            end
            
            Q_bar = mean(Fs,3);
            Q_norm = zeros(1,N_subs);
            
            % find the Q function nearest to the mean
            for i = 1:N_subs
                Q_norm(i) = sum((Fs(:,:,i) - Q_bar).^2 .* obj.A, 'all');
                disp(i)
            end
            
            idx = find(Q_norm == min(Q_norm));
            Q_mu = Fs(:,:,idx);
            
            % iterate closer the the karcher median
            for iter = 1:iters
                vv = zeros(size(Fs));
                distance = zeros(N_subs,1);
                    
                for i = 1:N_subs
                    tmpQ = Fs(:,:,i);
                    
                    tmpTheta = sum(tmpQ(:) .* Q_mu(:) .* obj.A(:));
                    
                    if 1 - abs(tmpTheta) < 1e-14
                        tmpTheta = sign(tmpTheta);
                    end

                    distance(i) = acos(tmpTheta);
                    
                    if distance(i) > 0
                        vv(:,:,i) = ((distance(i) / sin(distance(i))) * (tmpQ - cos(distance(i))*Q_mu)) / distance(i);
                    end
                end
                
                v_bar = sum(vv,3) / sum(1./distance(distance > 0));
                tmp = sqrt(sum(v_bar(:) .* v_bar(:) .* obj.A(:)));
                Q_mu = (cos(0.2*tmp) * Q_mu) + (sin(0.2*tmp) * (v_bar / tmp));
                Q_mu = Q_mu / sqrt(sum(Q_mu(:) .* Q_mu(:) .* obj.A(:)));

                fprintf('Template norm: %0.4f\n', tmp);

                if tmp < 0.005
                    break
                end
            end
            
            template = Q_mu;
        end

        function [result,lh_warp,rh_warp,cost] = register(obj,F1,F2,kernel,varargin)            
            p = inputParser;
            addParameter(p, 'is_template', false, @islogical);
            addParameter(p, 'verbose', 0, @isnumeric);

            % parse optional variables
            parse(p, varargin{:});
            params = p.Results;

            F2 = copy(F2);
            lh_warp = SphericalWarp(obj.lh_grid,1e-10);
            rh_warp = SphericalWarp(obj.rh_grid,1e-10);      

            % initial functions      
            if params.is_template
                Q1 = F1;
            else                
                Q1 = sqrt(F1.evaluate(kernel, true));
            end
            
            Q2 = sqrt(F2.evaluate(kernel, true));
            
            % initial cost
            moving_img = Q2;
            FmM = (Q1 - moving_img);

            cost = zeros(obj.max_iters+1,1);
            cost(1) = sum(FmM(:).^2 .* obj.A(:));            

            P = length(lh_warp.V);
            P2 = length(lh_warp.V) + length(rh_warp.V);

            last_lh_warp = lh_warp;
            last_rh_warp = rh_warp;
            
            fprintf('Initial cost for Sub: %0.6f\n', cost(1))

            idx_a = 1:P;
            idx_b = (P+1):P2;

            % start registering
            for iter = 1:obj.max_iters  
                % calculate the derivative
                [dQ2e1, dQ2e2] = obj.interpolator.get_derivative(moving_img);
            
                % evaluate derivative of the cost function
                FmM = FmM .* obj.A;
            
                % ------------------------------------------------
                % compute the gradient for the LH warp
                a = 2 * sum(FmM(idx_a,idx_a) .* dQ2e1(idx_a,idx_a),2) + ...
                           sum(FmM(idx_a,idx_b) .* dQ2e1(idx_a,idx_b),2);

                b = 2 * sum(FmM(idx_a,idx_a) .* dQ2e2(idx_a,idx_a),2) + ...
                           sum(FmM(idx_a,idx_b) .* dQ2e2(idx_a,idx_b),2);

                c = sum(FmM(idx_a,:) .* moving_img(idx_a,:),2);

                lh_dH = -2 * (sum(a.*obj.lh_grid.basis(:,:,1)) + ...
                              sum(b.*obj.lh_grid.basis(:,:,2)) + ...
                              sum(c.*obj.lh_grid.laplacian));
            
                % calculate displacement in each basis direction
                lh_step_size = obj.delta / (norm(lh_dH) + 1e-15);
                lh_gamma = squeeze(sum(lh_dH .* obj.lh_grid.basis,2)); 
                
                % ------------------------------------------------
                a = 2 * sum(FmM(idx_b,idx_b) .* dQ2e1(idx_b,idx_b),2) + ...
                           sum(FmM(idx_b,idx_a) .* dQ2e1(idx_b,idx_a),2);

                b = 2 * sum(FmM(idx_b,idx_b) .* dQ2e2(idx_b,idx_b),2) + ...
                           sum(FmM(idx_b,idx_a) .* dQ2e2(idx_b,idx_a),2);

                c = sum(FmM(idx_b,:) .* moving_img(idx_b,:),2);

                rh_dH = -2 * (sum(a.*obj.rh_grid.basis(:,:,1)) + ...
                              sum(b.*obj.rh_grid.basis(:,:,2)) + ...
                              sum(c.*obj.rh_grid.laplacian));           
            
                % calculate displacement in each basis direction
                rh_step_size = obj.delta / (norm(rh_dH) + 1e-15);
                rh_gamma = squeeze(sum(rh_dH .* obj.rh_grid.basis,2));                 

                % ------------------------------------------------   
                % compose the new warps with all previous warps                
                lh_warp = lh_warp.compose_warp(lh_step_size .* lh_gamma);                 
                rh_warp = rh_warp.compose_warp(rh_step_size .* rh_gamma);
                
                F2.warp_connectome(lh_warp, rh_warp);                
                moving_img = sqrt(F2.evaluate(kernel, true));
                
                % evaluate the new cost
                FmM = Q1 - moving_img; 
                cost(iter+1) = sum(FmM(:).^2 .* obj.A(:));
                
                if ((cost(iter) - cost(iter+1)) < obj.threshold)       
                   lh_warp = last_lh_warp;
                   rh_warp = last_rh_warp;
                   fprintf('Converged (increased cost) %d: %0.6f -> %0.6f\n', iter, cost(1), cost(iter))            
                   break
                end
            
                last_lh_warp = lh_warp;
                last_rh_warp = rh_warp;                           
            
                % show progress if asked
                if (params.verbose > 0)
                    fig = figure(params.verbose);
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
            end                      
           
            % evaluate function after warping
            F2.warp_connectome(lh_warp, rh_warp);                
            result = F2;

            cost = cost(1:iter+1);
        end
    end
end
