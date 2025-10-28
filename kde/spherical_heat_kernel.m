function kernel = spherical_heat_kernel(lh_grid,rh_grid,sigma,H)
    % calculate cutoff distance (ensures no negative connections)
    epsilon = 0.001;
    cross_harmonics = zeros(H,10000);
    
    Xs = linspace(1,-1,10000);
    for h = 1:H
        Ph = legendre(h,Xs,'norm');
        cross_harmonics(h,:) = Ph(1,:);
    end
    
    Hs = 1:H;
    for x = 1:10000    
        tmp = exp(-sigma * (Hs .* (Hs+1) + Hs' .* (Hs' + 1))) .* ((2*Hs)+1) .* ((2*Hs')+1) .* (cross_harmonics(:,1) * cross_harmonics(:,x)');
        tmp = sum(tmp(:));
    
        if tmp < epsilon       
            cutoff_dst = Xs(x-1);
            disp(cutoff_dst)
            break
        end
    end

    lh_kernel = calc_kernel(lh_grid,sigma,H,cutoff_dst);
    rh_kernel = calc_kernel(rh_grid,sigma,H,cutoff_dst);

    kernel = [lh_kernel, zeros(size(lh_kernel,1),size(rh_kernel,2));
              zeros(size(rh_kernel,1),size(lh_kernel,2)), rh_kernel];

    % normalise so rows add up to one
    kernel = sparse(kernel ./ sum(kernel,2).');
end

function kernel = calc_kernel(grid,sigma,H,cutoff_dst)
    % angles between all vertices
    X = max(min(grid.V * grid.V.',1),-1);
    
    % calculate the kernel for the grid
    kernel = zeros(size(X));
    
    for h = 1:H
        Ph = legendre(h,X,'norm'); 
        kernel = kernel + ((2*h) + 1) * exp(-h*(h+1)*sigma) .* squeeze(Ph(1,:,:));
        disp(h)
    end

    kernel(X < cutoff_dst) = 0;
end