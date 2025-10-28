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

        function template = get_template(~, Fs, kernel, iters)
            N_subs = size(Fs,1);
                        
            % find the mean Q function            
            for i = 1:N_subs
                tmp = sqrt(Fs{i}.evaluate(kernel,true));                

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
                Q_norm(i) = sum((sqrt(Fs{i}.evaluate(kernel,true)) - Q_bar).^2, 'all');
                disp(i)
            end
            
            idx = find(Q_norm == min(Q_norm),1);
            Q_mu = sqrt(Fs{idx}.evaluate(kernel,true));
            
            % iterate closer the the karcher median
            for iter = 1:iters
                vv = zeros(size(Q_mu,1),size(Q_mu,1),size(Fs,1));
                distance = zeros(N_subs,1);
                    
                for i = 1:N_subs
                    tmpQ = sqrt(Fs{i}.evaluate(kernel,true));
                    
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
           
            P = length(obj.lh_grid.V);
            P2 = length(obj.lh_grid.V) + length(obj.lh_grid.V);

            idx_a = 1:P;
            idx_b = (P+1):P2;
            
            lh_vt = 0;            
            rh_vt = 0;      

            % initial cost
            FmM = (Q1 - Q2);
            cost = zeros(obj.max_iters+1,1);
            cost(1) = sum(FmM(:).^2);            

            fprintf('Initial cost for Sub: %0.6f\n', cost(1))

            % start registering
            for iter = 1:obj.max_iters  
                % calculate the derivative
                [dQ2e1, dQ2e2] = obj.interpolator.get_derivative(Q2);            
            
                % ------------------------------------------------
                % compute the gradient for the LH warp
                % ------------------------------------------------
                a = 2*sum(FmM(idx_a,[idx_a,idx_b]) .* dQ2e1(idx_a,[idx_a,idx_b]),2) .* obj.lh_grid.basis(:,:,1);               
                b = 2*sum(FmM(idx_a,[idx_a,idx_b]) .* dQ2e2(idx_a,[idx_a,idx_b]),2) .* obj.lh_grid.basis(:,:,2);
                c = sum(FmM(idx_a,[idx_a,idx_b]) .* Q2(idx_a,[idx_a,idx_b]),2) .* obj.lh_grid.laplacian;
               
                lh_dH = 2 * sum(a+b+c,1);  
                lh_dH = squeeze(sum(lh_dH .* obj.lh_grid.basis,2));               
                
                % calculate displacement in each basis direction (using
                % momentum gradient descent to avoid local minimums)

                % TODO: CONSIDER CONSTRAINING MAX STEP SIZE TO KEEP
                % DIFFEOMORPHISMS SMOOTH (NO TRIANGLE FOLDS)
                lh_vt = 0.9*lh_vt + 0.1*lh_dH;                
                lh_step_size = -obj.delta .* lh_vt;                
                
                % ------------------------------------------------
                % compute the gradient for the RH warp
                % ------------------------------------------------
                a = 2*sum(FmM(idx_b,[idx_b,idx_a]) .* dQ2e1(idx_b,[idx_b,idx_a]),2) .* obj.rh_grid.basis(:,:,1);               
                b = 2*sum(FmM(idx_b,[idx_b,idx_a]) .* dQ2e2(idx_b,[idx_b,idx_a]),2) .* obj.rh_grid.basis(:,:,2);
                c = sum(FmM(idx_b,[idx_b,idx_a]) .* Q2(idx_b,[idx_b,idx_a]),2) .* obj.rh_grid.laplacian;
               
                rh_dH = 2 * sum(a+b+c,1);  
                rh_dH = squeeze(sum(rh_dH .* obj.rh_grid.basis,2));                

                % calculate displacement in each basis direction (using
                % momentum gradient descent to avoid local minimums)

                % TODO: CONSIDER CONSTRAINING MAX STEP SIZE TO KEEP
                % DIFFEOMORPHISMS SMOOTH (NO TRIANGLE FOLDS)
                rh_vt = 0.9*rh_vt + 0.1*rh_dH;                
                rh_step_size = -obj.delta .* rh_vt;                

                % ------------------------------------------------   
                % compose the new warps with all previous warps                
                % ------------------------------------------------   
                lh_warp.compose_warp(lh_step_size);                 
                rh_warp.compose_warp(rh_step_size);
                
                F2.warp_connectome(lh_warp, rh_warp);                
                Q2 = sqrt(F2.evaluate(kernel, true));
                
                % evaluate the new cost
                FmM = (Q1 - Q2); 
                cost(iter+1) = sum(FmM(:).^2);
                
                % TODO: CONSIDER LITTLE CHANGE IN COST FUNCTION OVER
                % MULTIPLE ITERATIONS TO BE THE CUTOFF TO ALLOW FOR
                % MOMENTUM OVERRUN BUT STILL HAVE A RELATIVLY LARGE
                % THRESHOLD (i.e. 1e-3 instead of 1e-5)
                if (abs(cost(iter) - cost(iter+1)) < obj.threshold)       
                   fprintf('Converged (increased cost) %d: %0.6f -> %0.6f\n', iter, cost(1), cost(iter))            
                   break
                end                                  
            
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
