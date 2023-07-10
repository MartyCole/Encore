function [dYlm,ddYlm] = Ylm_2nd_derivative(l, theta, phi)

% Make sure that if x is a 
% vector, it is a row vector
if size(theta,2) == 1
    theta = theta';
end
if size(phi,2) == 1
    phi = phi';
end

[ddPlmdx, dPlmdx, Plm] = legendre_2nd_derivative(l,theta);

m = (0:l)';
theta = permute(theta, [3 1:2]);
phi = permute(phi, [3 1:2]);
    
sphi = (sin(m .* phi));
cphi = (cos(m .* phi));
stheta = (max(sin(theta), 0.00001));

dYlm = zeros([2,size(dPlmdx)]);
ddYlm = zeros([4,size(ddPlmdx)]);

% Place results in tensor Re(Y_0), Re(Y_1), Im(Y_1), Re(Y_2), Im(Y_2), ...
p = 2*(l+1)-1;

dYlm(1,[1,2:2:p],:,:) = dPlmdx .* cphi;
dYlm(1,3:2:p,:,:) = dPlmdx(2:end,:,:) .* sphi(2:end,:,:);
dYlm(2,[1,2:2:p],:,:) = Plm .* -(m .* sphi) ./ stheta;
dYlm(2,3:2:p,:,:) = Plm(2:end,:,:) .* (m(2:end) .* cphi(2:end,:,:)) ./ stheta;

ddYlm(1,[1,2:2:p],:,:) = -l*(l + 1) * Plm .* cphi;
ddYlm(1,3:2:p,:,:) = -l*(l + 1) * Plm(2:end,:,:) .* sphi(2:end,:,:);
end