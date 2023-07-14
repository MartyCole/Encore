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

        function [result,lh_warp,rh_warp,cost] = register(obj,F1,F2) 
            lh_warp = SphericalWarp(obj.lh_grid,1e-10);
            rh_warp = SphericalWarp(obj.rh_grid,1e-10);

            % initial functions            
            Q1 = sqrt(F1 / sum(F1(:) .* obj.A(:)));
            Q2 = sqrt(F2 / sum(F2(:) .* obj.A(:)));
            
            % initial cost
            moving_img = Q2;
            FmM = (Q1 - moving_img);
            last_cost = sum(FmM(:).^2 .* obj.A(:));

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
            
                a = sum(FmM .* (2*dQ2e1), 1);
                b = sum(FmM .* (2*dQ2e2), 1);
                c = sum(FmM .* moving_img, 1);
                
                % ------------------------------------------------
                % compute the gradient for the LH warp
                lh_dH = 2 * (a(1:P) * obj.lh_grid.basis(:,:,1) + ...
                             b(1:P) * obj.lh_grid.basis(:,:,2) + ...
                             c(1:P) * obj.lh_grid.laplacian);               
            
                % calculate displacement in each basis direction
                lh_gamma = squeeze(sum(lh_dH .* obj.lh_grid.basis,2));            
                lh_step_size = min(obj.delta / max(vecnorm(lh_gamma')), 1);

                % ------------------------------------------------
                % compute the gradient for the RH warp
                rh_dH = 2 * (a((P+1):end) * obj.rh_grid.basis(:,:,1) + ...
                             b((P+1):end) * obj.rh_grid.basis(:,:,2) + ...
                             c((P+1):end) * obj.rh_grid.laplacian);               
            
                % calculate displacement in each basis direction
                rh_gamma = squeeze(sum(rh_dH .* obj.rh_grid.basis,2)); 
                rh_step_size = min(obj.delta / max(vecnorm(rh_gamma')), 1);    

                % ------------------------------------------------   
                % compose the new warps with all previous warps
                step_size = min(lh_step_size,rh_step_size);
                
                lh_warp = lh_warp.compose_warp(step_size .* lh_gamma);                 
                rh_warp = rh_warp.compose_warp(step_size .* rh_gamma);

                % evaluate function after warping
                moving_img = obj.concon.evaluate_Q(Q2,lh_warp,rh_warp);
                               
                % evaluate the new cost
                FmM = Q1 - moving_img; 
                cost = sum(FmM(:).^2 .* obj.A(:));
               
                if ((last_cost - cost) < obj.threshold)       
                   lh_warp = last_lh_warp;
                   rh_warp = last_rh_warp;
                   fprintf('Converged (increased cost) %d: %0.6f\n', iter, last_cost)            
                   break
                end
            
                last_lh_warp = lh_warp;
                last_rh_warp = rh_warp;
                last_cost = cost;                 
            
                % print progress
                if (mod(iter,10) == 0)       
                    fprintf('Iteration %d cost: %0.6f\n', iter, cost); 
                end
            end 
            
            result = obj.concon.evaluate(F2,lh_warp,rh_warp);            
        end
    end
end