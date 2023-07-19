function [b,laplacian,theta,phi,psi_norm] = tangent_basis(l, theta, phi, dA)

%% Generate psi functions
psi = zeros(size(theta,1),size(theta,2),(l+1)^2 - 1,2);
dpsi = zeros(size(theta,1),size(theta,2),(l+1)^2 - 1,4);
idx = 1;
for i=1:l
    len = 2*(i+1) - 1;
    % Note: MATLAB computes faster if evaluating over lower indices
    [dYlm,ddYlm] = Ylm_2nd_derivative(i,theta,phi);
    psi(:,:,idx:(idx+len-1),:) = permute(dYlm,[3,4,2,1]);
    dpsi(:,:,idx:(idx+len-1),:) = permute(ddYlm,[3,4,2,1]);
    idx = idx + len;
end

%% Find normalisation constants

% Perform double integral in spherical coordinates 
psi_norm = zeros(size(psi,3),1);
for i = 1:size(psi,3)
     psi_norm(i) = sqrt(sum((psi(:,:,i,1).^2 + psi(:,:,i,2).^2) .* dA, 'all'));         
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