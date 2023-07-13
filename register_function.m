function [result,warp] = register_function(grid,F1,F2,delta,iters)

% registration variables
interpolator = SphericalInterpolator(grid);
sph_derivative = SphericalDerivative(grid, 1e-10, 'tangent');

% initial warp    
warp = SphericalWarp(grid);
last_warp = warp;

% cached data
A = grid.A * grid.A';

% initial functions
Q1 = sqrt(F1 / sum(F1(:) .* A(:)));
Q2 = sqrt(F2 / sum(F2(:) .* A(:)));

% initial cost
Q_reg = Q2;
FmM = (Q1 - Q_reg);
last_cost = sum(FmM(:).^2 .* A(:));

fprintf('Initial cost for Sub: %0.6f\n', last_cost)

% start registering
for iter = 1:iters  
    % calculate the derivative
    [dQ2e1, dQ2e2] = sph_derivative.get_derivative(Q_reg);

    % evaluate derivative of the cost function
    FmM = FmM .* A;

    a = sum(FmM .* (2*dQ2e1), 1);
    b = sum(FmM .* (2*dQ2e2), 1);
    c = sum(FmM .* Q_reg, 1);
      
    dH = 2 * ((a * grid.basis(:,:,1)) + (b * grid.basis(:,:,2)) + (c * grid.laplacian));

    if norm(dH,'fro') < 1e-6
        fprintf('Converged (small norm) %d: %0.6f\n', iter, cost)            
        break   
    end    

    % calculate displacement in each basis direction
    gamma = squeeze(sum(dH .* grid.basis,2));            
    step_size = min(0.02 / max(vecnorm(gamma')), delta);
    step_size = delta;

    % compose the new warp with all previous warps
    warp = sph_derivative.compose_warp(step_size .* gamma, warp);

    % evaluate function after warping
    interpolator.update_query_points(warp.V);
    Q_reg = interpolator.evaluate_2d(Q2) .* (sqrt(warp.J) * sqrt(warp.J).');
    Q_reg = (Q_reg + Q_reg.') / 2;
    
    Q_norm = sum(Q_reg(:).^2 .* A(:));
    Q_reg = Q_reg / sqrt(Q_norm);
    
    % evaluate the new cost
    FmM = Q1 - Q_reg; 
    cost = sum(FmM(:).^2 .* A(:));
   
    if ((last_cost - cost) < 0)       
       warp = last_warp;
       fprintf('Converged (increased cost) %d: %0.6f\n', iter, last_cost)            
       break
    end

    last_warp = warp;
    last_cost = cost; 

    % print progress
    if (mod(iter,10) == 0)       
        fprintf('Iteration %d cost: %0.6f, with norm: %0.4f\n', iter, cost, Q_norm); 
    end
end 

result = interpolator.evaluate_2d(F2) .* (warp.J * warp.J.');

end