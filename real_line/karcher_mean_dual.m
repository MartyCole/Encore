function template = karcher_mean_dual(Q,iters,v,dv)
% initialisation
T = size(Q,1) / 2;
N = size(Q,3);
warp_hat_1 = zeros(T,N);
warp_hat_2 = zeros(T,N);

% find the mean Q function
Q_bar = mean(Q,3);
Q_norm = zeros(1,N);

% find the Q function nearest to the mean
for i = 1:N
    Q_norm(i) = norm(Q(:,:,i) - Q_bar, 'fro');
    warp_hat_1(:,i) = linspace(0,1,T);
    warp_hat_2(:,i) = linspace(0,1,T);
end

idx = (Q_norm == min(Q_norm));
Q_mu = Q(:,:,idx);

% iterate closer the the karcher mean
for iter = 1:iters
    vv = zeros(2*T,2*T,N);
    for i = 1:N
        tmpQ = Q(:,:,i);
        
        theta = acos(trapz(trapz((tmpQ.*Q_mu))));
        
        if theta > 0
            vv(:,:,i) = (theta / sin(theta)) * (tmpQ - cos(theta)*Q_mu);
        end
    end
    
    v_bar = mean(vv,3);    
    tmp = sqrt(trapz(trapz(v_bar.^2)));
    Q_mu = (cos(0.2*tmp) * Q_mu) + (sin(0.2*tmp) * (v_bar / tmp));
    Q_mu = Q_mu / norm(Q_mu,'fro');
end

% % find the mean warp
% for i = 1:N
%    [warp_hat_1(:,i),warp_hat_2(:,i),~] = grad_search(Q_mu,Q(:,:,i),...
%                                          v,dv,iters,1,0.1); 
% end
% 
% phi_1 = zeros(T,N);
% 
% for i = 1:N
%     phi_1(:,i) = sqrt(gradient(warp_hat_1(:,i), 1 / (T-1)));
% end
% 
% phi_mu = phi_1(:,1);
% vec = zeros(T,N);
% 
% for iter = 1:100
%     for i = 1:N
%         theta = trapz(phi_mu.*phi_1(:,i))*(1/(T-1));
%         theta = max(0,min(1,theta));
% 
%         theta = acos(theta);
% 
%         if theta > 0
%             vec(:,i) = (theta/sin(theta))*(phi_1(:,i)-cos(theta)*phi_mu);
%         else
%             vec(:,i) = 0;
%         end
%     end
% 
%     v_bar = mean(vec,2);
%     norm_v = norm(v_bar);
% 
%     if norm_v < 0.0001
%         break;
%     end
%     
%     phi_mu = cos(0.1*norm_v)*phi_mu + sin(0.1*norm_v)*(v_bar/norm_v);
% end
% 
% warp_bar = cumtrapz(phi_mu.^2);
% warp_bar = warp_bar / warp_bar(end);
% 
% warp_bar_inv = invert_gamma(warp_bar)';

% template = linear_interp2D_q(Q_mu, warp_bar_inv);
template = Q_mu;
end