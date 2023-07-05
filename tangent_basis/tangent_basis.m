function [b,laplacian,theta,phi,psi_norm] = tangent_basis(l, varargin)

% TODO: Function information
if nargin == 2
    n = varargin{1};
    
    % Generate a grid to evaluate over
    theta = repmat(linspace(0,pi,n),n,1);
    phi = repmat(linspace(0,2*pi,n)',1,n);
    
    % Stepsize to calculate norm
    dtheta = pi/(n-1);
    dphi = (2*pi)/(n-1);
elseif nargin == 4
    % Find grid to evaluate over
    theta = varargin{1};
    phi = varargin{2};
    psi_norm = varargin{3};
    
    if size(psi_norm,1) ~= ((l+1)^2 - 1)
        error('Normalising constant does not match number of psi functions');
    end
else
    error('Wrong number of arguments');
end

%% Generate psi functions
psi = zeros(size(theta,1),size(theta,2),(l+1)^2 - 1,2);
dpsi = zeros(size(theta,1),size(theta,2),(l+1)^2 - 1,4);
idx = 1;
for i=1:l
    len = 2*(i+1) - 1;
    % Note: MATLAB computes faster if evaluating over lower indices
    %[dYlm] = Ylm_derivative(i,theta,phi);
    [dYlm,ddYlm] = Ylm_2nd_derivative(i,theta,phi);
    psi(:,:,idx:(idx+len-1),:) = permute(dYlm,[3,4,2,1]);
    dpsi(:,:,idx:(idx+len-1),:) = permute(ddYlm,[3,4,2,1]);
    idx = idx + len;
end

%% Find normalisation constants

% Perform double integral in spherical coordinates 
% TODO: Try to derive the constant to skip this step
if nargin == 2
    tol = 1e-12;
    sx = sin(linspace(0,pi,n)');
    sx(sx <= tol) = 0;

    psi_norm = zeros(size(psi,3),1);

    for i = 1:size(psi,3)
         theta_norm = trapz(trapz(psi(:,:,i,1).*psi(:,:,i,1) .* sx)) * dtheta * dphi;
         phi_norm = trapz(trapz(psi(:,:,i,2).*psi(:,:,i,2) .* sx)) * dtheta * dphi;

         psi_norm(i) = sqrt(theta_norm + phi_norm);
    end
end

% Apply the normalisation constants
for i = 1:size(psi,3)
     if psi_norm(i) > 0
        psi(:,:,i,1) = psi(:,:,i,1) ./ psi_norm(i);
        psi(:,:,i,2) = psi(:,:,i,2) ./ psi_norm(i);     
        dpsi(:,:,i,1) = dpsi(:,:,i,1) ./ psi_norm(i);
     end
end

%% Compute the complete basis system 
b = zeros(size(theta,1),size(theta,2),2*(l+1)^2 - 2,2);
laplacian = zeros(size(theta,1),size(theta,2),2*(l+1)^2 - 2);

b(:,:,1:(idx-1),1:2) = psi(:,:,1:(idx-1),:);
b(:,:,idx:end,1) =  psi(:,:,1:(idx-1),2);
b(:,:,idx:end,2) = -psi(:,:,1:(idx-1),1);

laplacian(:,:,1:(idx-1)) = dpsi(:,:,1:(idx-1));
laplacian(:,:,idx:end) = 0;

end