function [gamma_1, gamma_2, Q_hat] = grad_search_dual(Q1,Q2,v,dv,iters,lambda,delta)

P = size(Q1,1) / 2;

LH_Q1 = Q1(1:P,1:P);
LH_Q2 = Q2(1:P,1:P);

RH_Q1 = Q1((P+1):end,(P+1):end);
RH_Q2 = Q2((P+1):end,(P+1):end);

BH_Q1 = Q1(1:P,(P+1):end);
BH_Q2 = Q2(1:P,(P+1):end);

old_LH_Q2 = LH_Q2;
old_RH_Q2 = RH_Q2;
old_BH_Q2 = BH_Q2;

gamma_1 = linspace(0,1,P)';
gamma_2 = linspace(0,1,P)';

last_cost = inf;

for iter = 1:iters
    % Iterate the left hemisphere
    [dLH,dLHRH] = new_cost(LH_Q1,LH_Q2,BH_Q1,BH_Q2,v,dv);
    [dRH,dRHLH] = new_cost(RH_Q1,RH_Q2,BH_Q1',BH_Q2',v,dv);
    %[dLHRH,~] = new_cost(BH_Q1,BH_Q2,BH_Q1,BH_Q2,v,dv);
    %[dRHLH,~] = new_cost(BH_Q1',BH_Q2',BH_Q1,BH_Q2',v,dv);

    %dL = ((1-lambda)*dH) + (lambda * dK);
    dL = dLH + (lambda*dLHRH);

    test_gamma_1 = gamma_1;
    test_gamma_2 = gamma_2;

    if norm(dL) > 0.02
        new_gamma_1 = generate_gamma(dL,gamma_1,delta);
    
        if any(gradient(new_gamma_1, 1 / (P-1)) < 0)        
            break
        end
        
        test_gamma_1 = new_gamma_1;
        
        
        %LH_Q2 = (LH_Q2 + LH_Q2') / 2;
        %LH_Q2 = LH_Q2 / norm(LH_Q2);
    
        %BH_Q2 = linear_interp2D_q(old_BH_Q2, gamma_1, gamma_2);
        %BH_Q2 = BH_Q2 / norm(BH_Q2);
    end

    % Iterate the right hemisphere
    %[dH,dK] = new_cost(RH_Q1,RH_Q2,BH_Q1',BH_Q2',v,dv);

    %dL = ((1-lambda)*dH) + (lambda * dK);
    dL = dRH + (lambda*dLHRH);
    
    if norm(dL) > 0.02
        new_gamma_2 = generate_gamma(dL,gamma_2,delta);
    
        if any(gradient(new_gamma_2, 1 / (P-1)) < 0)        
            break
        end
        
        test_gamma_2 = new_gamma_2;           
    end

    test_LH_Q2 = linear_interp2D_q(old_LH_Q2, test_gamma_1, test_gamma_1);
    test_RH_Q2 = linear_interp2D_q(old_RH_Q2, test_gamma_2, test_gamma_2);
    test_BH_Q2 = linear_interp2D_q(old_BH_Q2, test_gamma_1, test_gamma_2);
    %BH_Q2 = BH_Q2 / norm(BH_Q2);           

    test_Q = [test_LH_Q2,test_BH_Q2;test_BH_Q2',test_RH_Q2];
    test_Q = (test_Q + test_Q') / 2; 
    test_Q = test_Q / norm(test_Q); 
    
    cost = sum(sum((Q1 - test_Q).^2)); 

    if cost < last_cost
        Q_hat = test_Q;
        LH_Q2 = Q_hat(1:P,1:P);      
        RH_Q2 = Q_hat((P+1):end,(P+1):end);
        BH_Q2 = Q_hat(1:P,(P+1):end);  

        gamma_1 = test_gamma_1;
        gamma_2 = test_gamma_2;

        last_cost = cost;
    else
        break;
    end   
end

Q_hat = [LH_Q2,BH_Q2;BH_Q2',RH_Q2];

end

function [dH, dK] = new_cost(Q1, Q2, QI1, QI2, v, dv)
    P = size(Q1,1);
    R = size(v,2);
        
    % Cost for intra-connections
    dH = zeros(P,1);

    [dQ2x,dQ2y] = gradient(Q2, 1 / (P-1));     
    Q1mQ2 = Q1 - Q2;
    Q1Q2 = (Q1mQ2) .* Q2;

    for i = 1:R
        dQv = sum(v(:,i)' * ((Q1mQ2) .* dQ2x)) + sum(((Q1mQ2) .* dQ2y) * v(:,i));
        
        Qvx = 0.5 * sum(dv(:,i)' * Q1Q2);
        Qvy = 0.5 * sum(Q1Q2 * dv(:,i));
        
        %dH(:) = 2 * sum((Q1mQ2 .* (dQ2x' + dQ2y)) * v(:,i) + (Q1Q2 * dv(:,i))) * v(:,i);

        dH(:) = dH(:) + 2 * (dQv + Qvx + Qvy) * v(:,i); 
    end

    % Cost for inter-connections
    dK = zeros(P,1);
    
    [dQ2x,~] = gradient(QI2, 1 / (P-1));     
    Q1mQ2 = QI1 - QI2;
    Q1Q2 = (Q1mQ2) .* QI2;

    for i = 1:R
        dQv = sum(v(:,i)' * ((Q1mQ2) .* dQ2x));        
        Qvx = 0.5 * sum(dv(:,i)' * Q1Q2); 
        
        dK(:) = dK(:) + 2 * (dQv + Qvx) * v(:,i); 
    end
end

function new_gamma = generate_gamma(dH,gamma,eps)
    P = length(gamma);
     
    while eps > 1e-6
        test_warp = linspace(0,1,P)' + eps*dH;

        if any(gradient(test_warp, 1 / (P-1)) < 0)
            eps = eps / 2;
        else
            break;
        end
    end
    
    if eps <= 1e-6
        error('Cannot find a small enough epsilon');       
    end

    new_gamma = interp1(linspace(0,1,P),gamma,test_warp);   
    new_gamma(1) = 0;    
    new_gamma(end) = 1;
end