function [start_points,end_points] = random_connectome(N,kappa1,mu1,kappa2,mu2)

mu1_cart = sphere_to_cart(mu1(1),mu1(2));
mu2_cart = sphere_to_cart(mu2(1),mu2(2));

% Generate points from which connections start and end
lr_sp = rand_von_mises(kappa1, 3, mu1_cart, N);
lr_ep = rand_von_mises(kappa1, 3, mu2_cart, N);
lh_sp = rand_von_mises(kappa2, 3, mu1_cart, N);
lh_ep = rand_von_mises(kappa2, 3, mu1_cart, N);
rh_sp = rand_von_mises(kappa2, 3, mu2_cart, N);
rh_ep = rand_von_mises(kappa2, 3, mu2_cart, N);

% Mask LH points to be only on LH
mask = (lh_sp(1,:) < 0) | (lh_ep(1,:) < 0);
lh_sp = lh_sp(:,~mask);
lh_ep = lh_ep(:,~mask);

% Mask RH points to be only on RH
mask = (rh_sp(1,:) > 0) | (rh_ep(1,:) > 0);
rh_sp = rh_sp(:,~mask);
rh_ep = rh_ep(:,~mask);

% Concatenate all points
start_points = [lr_sp,lh_sp,rh_sp];
end_points = [lr_ep,lh_ep,rh_ep];

end

