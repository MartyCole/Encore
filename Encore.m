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
        concon
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
             
            obj.concon = Concon(obj.lh_grid, obj.rh_grid, 1e-10);

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

        function [result,lh_warp,rh_warp,cost] = register(obj,F1,F2,varargin)
            p = inputParser;
            addParameter(p, 'is_template', false, @islogical);
            addParameter(p, 'verbose', 0, @isnumeric);

            % parse optional variables
            parse(p, varargin{:});
            params = p.Results;

            lh_warp = SphericalWarp(obj.lh_grid,1e-10);
            rh_warp = SphericalWarp(obj.rh_grid,1e-10);

            % initial functions      
            if params.is_template
                Q1 = F1;
            else                
                Q1 = sqrt(F1 / sum(F1(:) .* obj.A(:)));
            end

            Q2 = sqrt(F2 / sum(F2(:) .* obj.A(:)));
            
            % initial cost
            moving_img = Q2;
            FmM = (Q1 - moving_img);
            last_cost = sum(FmM(:).^2 .* obj.A(:));
            init_cost = last_cost;

            P = length(lh_warp.V);
            last_lh_warp = lh_warp;
            last_rh_warp = rh_warp;
            
            fprintf('Initial cost for Sub: %0.6f\n', last_cost)

            % start registering
            for iter = 1:obj.max_iters  
                % calculate the derivative
                [dQ2e1, dQ2e2] = obj.concon.get_derivative(moving_img);
            
                % evaluate derivative of the cost function
                FmM = FmM .* obj.A;
            
                a = sum(FmM .* (2*dQ2e1), 2).';
                b = sum(FmM .* (2*dQ2e2), 2).';
                c = FmM .* moving_img;
                
                % ------------------------------------------------
                % compute the gradient for the LH warp
                lh_dH = 2 * (a(1:P) * obj.lh_grid.basis(:,:,1) + ...
                             b(1:P) * obj.lh_grid.basis(:,:,2) + ...
                             sum(c(1:P,:),2).' * obj.lh_grid.laplacian);               
            
                % calculate displacement in each basis direction
                lh_step_size = obj.delta / (norm(lh_dH) + 1e-15);
                lh_gamma = squeeze(sum(lh_dH .* obj.lh_grid.basis,2)); 
                
                % ------------------------------------------------
                % compute the gradient for the RH warp
                rh_dH = 2 * (a((P+1):end) * obj.rh_grid.basis(:,:,1) + ...
                             b((P+1):end) * obj.rh_grid.basis(:,:,2) + ...
                             sum(c((P+1):end,:),2).' * obj.rh_grid.laplacian);               
            
                % calculate displacement in each basis direction
                rh_step_size = obj.delta / (norm(rh_dH) + 1e-15);
                rh_gamma = squeeze(sum(rh_dH .* obj.rh_grid.basis,2));                 

                % ------------------------------------------------   
                % compose the new warps with all previous warps                
                lh_warp = lh_warp.compose_warp(lh_step_size .* lh_gamma);                 
                rh_warp = rh_warp.compose_warp(rh_step_size .* rh_gamma);

                % evaluate function after warping
                moving_img = obj.concon.evaluate_Q(Q2,lh_warp,rh_warp);
                               
                % evaluate the new cost
                FmM = Q1 - moving_img; 
                cost = sum(FmM(:).^2 .* obj.A(:));
                
                if ((last_cost - cost) < obj.threshold)       
                   lh_warp = last_lh_warp;
                   rh_warp = last_rh_warp;
                   fprintf('Converged (increased cost) %d: %0.6f -> %0.6f\n', iter, init_cost, last_cost)            
                   break
                end
            
                last_lh_warp = lh_warp;
                last_rh_warp = rh_warp;
                last_cost = cost;                 
            
                % print progress
                if (mod(iter,10) == 0)       
                    fprintf('Iteration %d cost: %0.6f\n', iter, cost);                     
                end

                if (params.verbose > 0)
                    figure(params.verbose)
                    lh_warp.plot('Estimated LH-Warp')
                    clim([0 2])
                    view([0 90 0]);
                end
            end             

            result = obj.concon.evaluate(F2,lh_warp,rh_warp);              
        end
    end
end
