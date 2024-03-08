% The purpose of this script:
%   * Simulate registration of a population using gradient descent

addpath('../real_line')
clear variables;
close all;

%% Setup
% For reprodicibility
rng(02092022);

ITERS = 200;
DELTA = 0.1;

% Smooth connectivity functions before processing?
PRE_SMOOTH = true;
smooth_sigma = 5;

% Number of q functions to simulate
N = 80;

% Grid of size T over which to evaluate the functions
T = 200;
t = linspace(0,1,T);

% Number of tangent basis functions for warp functions
R = 120;

% Tolerance of small numbers
tol = 1e-10;

% Keeps track of figure ID
figID = 0;

% Default subject to plot
render_sample = 3;

%% Setup tangent basis functions
% Generate tangent basis functions, v, with rank R
r = 1:(R/2);
v = (1/(sqrt(2)*pi)) * [(1./r) .* sin(2*pi*r.*t'), ...
                        (1./r) .* (cos(2*pi*r.*t') - 1)];

v = [v(:,1:2:end),v(:,2:2:end)];
                  
% Generate the derivative of the basis functions
dv = (2/(sqrt(2))) * [ cos(2*pi*r.*t'), ...
                      -sin(2*pi*r.*t')];
                  
dv = [dv(:,1:2:end),dv(:,2:2:end)];

% set values near 0 to 0                    
v(abs(v) < tol) = 0;                    
dv(abs(dv) < tol) = 0;

%% Simulate 1D connectivity using gaussian endpoints
%    * Randomly generate endpoints for each subject using same parameters
%    * Add noise?
%    * Apply a random warp to each subject
%    * Generate a templat (Frechet mean)

% generate N random warps
orig_warps = random_noisy_warp1D(N, T, 1); 

M = zeros(T,T,N);
Q = zeros(T,T,N);
Q_unorm = zeros(T,T,N);

% Generate a population of N connectivity functions
for i = 1:N     
    [tmp_a, sp_a, ep_a] = simulate_connectivity1D(2500, T, 0.1,  0.9,  0.5, false);
    [tmp_b, sp_b, ep_b] = simulate_connectivity1D(2500, T, 0.1,  0.3,  0.5, false);
    [tmp_c, sp_c, ep_c] = simulate_connectivity1D(2500, T, 0.5,  0.6,  0.4, false);
    [tmp_d, sp_d, ep_d] = simulate_connectivity1D(2500, T, 0.2,  0.6,  0.2, false);
    [tmp_e, sp_e, ep_e] = simulate_connectivity1D(8000, T, 0.55, 0.95, 0.3, false);
    [tmp_f, sp_f, ep_f] = simulate_connectivity1D(2000, T, 0.4,  0.8,  0.6, false);
    
    sp = [sp_a;sp_b;sp_c;sp_d;sp_e;sp_f];
    ep = [ep_a;ep_b;ep_c;ep_d;ep_e;ep_f]; 

    % warp and smooth the matrices if specified
    M(:,:,i) = warp_connectivity1D(sp, ep, linspace(0,1,T), T);
    Q(:,:,i) = warp_connectivity1D(sp, ep, orig_warps(i,:)', T);
    
    if PRE_SMOOTH == true
        M(:,:,i) = imgaussfilt(M(:,:,i), smooth_sigma);
        Q(:,:,i) = imgaussfilt(Q(:,:,i), smooth_sigma);
    end
    
    % transform to the q space
    Q_unorm(:,:,i) = sqrt(M(:,:,i) / trapz(trapz(M(:,:,i))));
    Q(:,:,i) = sqrt(Q(:,:,i) / trapz(trapz(Q(:,:,i))));   
end

%% Plot warp functions
figID = figID + 1;
figure(figID)

plot(1:T, orig_warps(:,:)); hold on; plot(1:T, linspace(0,1,T)); hold off;
pbaspect([1,1,1]);
xlim([1, T])
xticklabels(linspace(0.1,1,10));
title('Original warp functions for population generation')

%% Plot warped subject
render_sample = 21;
figID = figID + 1;
f = figure(figID);

subplot(1,3,1);
pbaspect([1,1,1]);
imagesc(sqrt(M(:,:,render_sample)) / norm(sqrt(M(:,:,render_sample))));
pbaspect([1,1,1]);
xticklabels([]);
yticklabels([]);
title('Original connectivity')
subplot(1,3,2)
plot(1:T, orig_warps(render_sample, :), 'r-.')
pbaspect([1,1,1]);
refline(1/(T+1),0)
xlim([0, T])
xticks(linspace(0.1,1,10)*T);
xticklabels(linspace(0.1,1,10));
title('Warp')
subplot(1,3,3);
pbaspect([1,1,1]);
imagesc(Q(:,:,render_sample))
pbaspect([1,1,1]);
xticklabels([]);
yticklabels([]);
title('Warped connectivity')

f.Position = [50 500 1500 500];
sgtitle(sprintf('Warping connectivity for sample %i', render_sample))

%% Plot population warped connectivity
figID = figID + 1;
f = figure(figID);

imagesc(mean(M,3))
pbaspect([1,1,1]);
xticklabels([]);
yticklabels([]);
title('Element-wise mean connectivity after warping')

f.Position = [50 500 500 500];

%% Find Karcher mean of orbit Q

% initialisation
Q_hat = Q;
warp_hat = zeros(T,N);
last_norm = 0;

% find the mean Q function
Q_bar = mean(Q_hat,3);
Q_norm = zeros(1,N);

for i = 1:N
    Q_norm(i) = norm(Q_hat(:,:,i) - Q_bar, 'fro');
    warp_hat(:,i) = linspace(0,1,T);
end

idx = find(Q_norm == min(Q_norm));
Q_mu = Q(:,:,idx);

for iter = 1:200
    vv = zeros(200,200,N);
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

%% Find the mean warping function

for i = 1:N
   [~,warp_hat(:,i),~] = grad_search(Q_mu,Q_hat(:,:,i),...
                                         linspace(0,1,T),v,dv,ITERS,0.1); 
end

phi = zeros(T,N);

for i = 1:N
    phi(:,i) = sqrt(gradient(warp_hat(:,i), 1 / (T-1)));
end

phi_mu = phi(:,1);
vec = zeros(T,N);

for iter = 1:100
    for i = 1:N
        theta = trapz(phi_mu.*phi(:,i))*(1/(T-1));
        theta = max(0,min(1,theta));

        theta = acos(theta);

        if theta > 0
            vec(:,i) = (theta/sin(theta))*(phi(:,i)-cos(theta)*phi_mu);
        else
            vec(:,i) = 0;
        end
    end

    v_bar = mean(vec,2);
    norm_v = norm(v_bar);

    if norm_v < 0.00001
        break;
    end
    
    phi_mu = cos(0.1*norm_v)*phi_mu + sin(0.1*norm_v)*(v_bar/norm_v);
end

disp('done')

warp_bar = cumtrapz(phi_mu.^2);
warp_bar = warp_bar / warp_bar(end);

warp_bar_inv = invert_gamma(warp_bar)';

%% Define the template to register to

template = linear_interp2D_q(Q_mu, warp_bar_inv);

%% Register population to the template

Q_tilde = Q;
warp_tilde = zeros(T,N);

for i = 1:N
    [~,warp_tilde(:,i),Q_tilde(:,:,i)] = grad_search(template,Q(:,:,i),...
                                         linspace(0,1,T),v,dv,ITERS,DELTA);
    disp(i);
end

%% Nudge results back towards having identity warp for comparison

foo = zeros(T,T,N);

for i = 1:N
    foo(:,:,i) = sqrt(M(:,:,i));
    foo(:,:,i) = foo(:,:,i) / norm(foo(:,:,i));
end

[cost,mean_warp,~] = grad_search(mean(foo,3),mean(Q_hat,3),...
                                         linspace(0,1,T),v,dv,ITERS,DELTA);
                                     
Q_final = zeros(T,T,N);
fnl_warp = zeros(T,N);

for i = 1:N
   fnl_warp(:,i) = interp1(linspace(0,1,T),warp_hat(:,i),mean_warp);
   Q_final(:,:,i) = linear_interp2D_q(Q(:,:,i), fnl_warp(:,i)); 
   Q_final(:,:,i) = Q_final(:,:,i) / norm(Q_final(:,:,i));
end

%%
figID = figID + 1;

hb = 0.007;
for i = 1:N
f = figure(figID);
render_sample = i;
subplot(1,4,1)
plot(1:T, orig_warps(render_sample, :), 'r.', 'linewidth', 1) 
hold on
plot(1:T, interp1(fnl_warp(:,render_sample),linspace(0,1,T),linspace(0,1,T)), 'b-', 'linewidth', 1)
pbaspect([1,1,1]);
refline(1/(T+1),0)
xlim([0, T])
xticks(linspace(0.1,1,10)*T);
xticklabels(linspace(0.1,1,10));
hold off
title('Warp')
subplot(1,4,2);
pbaspect([1,1,1]);
imagesc(Q(:,:,render_sample))
clim([0,hb]);
pbaspect([1,1,1]);
xticks(linspace(0.1,1,10)*T);
xticklabels(linspace(0.1,1,10));
yticks(linspace(0.1,1,10)*T);
yticklabels(linspace(0.1,1,10));
set(gca,'YDir','normal')
title('Warped connectivity')
subplot(1,4,3);
pbaspect([1,1,1]);
imagesc(Q_tilde(:,:,render_sample))
clim([0,hb]);
pbaspect([1,1,1]);
xticks(linspace(0.1,1,10)*T);
xticklabels(linspace(0.1,1,10));
yticks(linspace(0.1,1,10)*T);
yticklabels(linspace(0.1,1,10));
set(gca,'YDir','normal')
title('Recovered connectivity')
subplot(1,4,4);
pbaspect([1,1,1]);
imagesc(foo(:,:,render_sample));
clim([0,hb]);
pbaspect([1,1,1]);
xticks(linspace(0.1,1,10)*T);
xticklabels(linspace(0.1,1,10));
yticks(linspace(0.1,1,10)*T);
yticklabels(linspace(0.1,1,10));
set(gca,'YDir','normal')
title('Original connectivity')

f.Position = [50 500 1500 500];
sgtitle(sprintf('Warping connectivity for sample %i', render_sample))
end
