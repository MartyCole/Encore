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

% Calculate normalisation term, C
m = (0:l)';
a1 = (2*l+1) / (4*pi);
a2 = factorial(l-m) ./ factorial(l+m);
C = sqrt(a1*a2);

theta = permute(theta, [3 1:2]);
phi = permute(phi, [3 1:2]);
    
sphi = (sin(m .* phi));
cphi = (cos(m .* phi));
stheta = (max(sin(theta), 0.00001));

dYlm = zeros([2,size(dPlmdx)]);
ddYlm = zeros([4,size(ddPlmdx)]);

% Place results in tensor Re(Y_0), Re(Y_1), Im(Y_1), Re(Y_2), Im(Y_2), ...
p = 2*(l+1)-1;

% TODO: Hopefully we can derive a normalisation term so that
%       we can skip the numerical calculation in the next step and we can
%       get function values for any arbitrary coordinates instead of having
%       to predefine a large grid and doing the double integration.
dYlm(1,[1,2:2:p],:,:) = C .* dPlmdx .* cphi;
dYlm(1,3:2:p,:,:) = C(2:end) .* dPlmdx(2:end,:,:) .* sphi(2:end,:,:);
dYlm(2,[1,2:2:p],:,:) = C .* Plm .* -(m .* sphi) ./ stheta;
dYlm(2,3:2:p,:,:) = C(2:end) .* Plm(2:end,:,:) .* (m(2:end) .* cphi(2:end,:,:)) ./ stheta;

% % % C .* dPlmdx .* cphi;
% % ddYlm(1,[1,2:2:p],:,:) = C .* ddPlmdx .* cphi; % dt/dt
% % ddYlm(3,[1,2:2:p],:,:) = C .* dPlmdx .* -(m .* sphi); % dt/tp
% % % C .* Plm .* -(m .* sphi);
% % ddYlm(2,[1,2:2:p],:,:) = C .* dPlmdx .* -(m .* sphi) ./ stheta.^2; % dp/dt
% % ddYlm(4,[1,2:2:p],:,:) = C .* Plm .* -(m.^2 .* cphi) ./ stheta.^2; % dp/dp
% % 
% % % C(2:end) .* dPlmdx(2:end,:,:) .* sphi(2:end,:,:);
% % ddYlm(1,3:2:p,:,:) = C(2:end) .* ddPlmdx(2:end,:,:) .* sphi(2:end,:,:); % dt/dt
% % ddYlm(3,3:2:p,:,:) = C(2:end) .* dPlmdx(2:end,:,:) .* (m(2:end) .* cphi(2:end,:,:)); %dt/tp
% % % C(2:end) .* Plm(2:end,:,:) .* (m(2:end) .* cphi(2:end,:,:));
% % ddYlm(2,3:2:p,:,:) = C(2:end) .* dPlmdx(2:end,:,:) .* (m(2:end) .* cphi(2:end,:,:)) ./ stheta.^2;
% % ddYlm(4,3:2:p,:,:) = C(2:end) .* Plm(2:end,:,:) .* (-m(2:end).^2 .* sphi(2:end,:,:)) ./ stheta.^2;

ddYlm(1,[1,2:2:p],:,:) = -l*(l + 1) * C .* Plm .* cphi;
ddYlm(1,3:2:p,:,:) = -l*(l + 1) * C(2:end,1) .* Plm(2:end,:,:) .* sphi(2:end,:,:);
end