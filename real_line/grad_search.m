function [cost, gamma, Q2] = grad_search(Q1,Q2,gamma,v,dv,iters,delta)

oldQ2 = Q2;
R = size(v,2);
T = length(gamma);

if ~all(gamma == linspace(0,1,T))
    Q2 = linear_interp2D_q(oldQ2, gamma);  
    Q2 = Q2 / norm(Q2);
end

for iter = 1:iters
    [dQ2x,dQ2y] = gradient(Q2, 1 / (T-1)); 
    
    dH = zeros(T,1);
    Q1mQ2 = Q1 - Q2;
    
    for i = 1:R
        dQv = sum(v(:,i)' * ((Q1mQ2) .* dQ2x)) + sum(((Q1mQ2) .* dQ2y) * v(:,i));

        Q1Q2 = (Q1mQ2) .* Q2;
        Qvx = 0.5 * sum(dv(:,i)' * Q1Q2);
        Qvy = 0.5 * sum(Q1Q2 * dv(:,i));
        
        dH(:) = dH(:) + 2 * (dQv + Qvx + Qvy) * v(:,i); 
    end

    eps = delta;
     
    while eps > 1e-6
        test_warp = linspace(0,1,T)' + eps*dH;

        if any(gradient(test_warp, 1 / (T-1)) < 0)
            eps = eps / 2;
        else
            break;
        end
    end
    
    if eps <= 1e-6
        error('Cannot find a small enough epsilon');       
    end
    
    old_gamma = gamma;
    gamma = interp1(linspace(0,1,T),gamma,test_warp);   
    gamma(1) = 0;    
    gamma(end) = 1;

    if any(gradient(gamma, 1 / (T-1)) < 0)
        gamma = old_gamma;
        break
    end
    
    Q2 = linear_interp2D_q(oldQ2, gamma);  
    Q2 = Q2 / norm(Q2);
    
    cost = sum(sum((Q1 - Q2).^2));
    
    %disp(cost)
    
    if norm(dH) < 0.0001
        disp('converged:')        
        break
    end
end
disp(cost)
end