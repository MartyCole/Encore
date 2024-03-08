% MATLAB R2018a
%
% FUNCTION NAME:
%   invert_gamma
%
% DESCRIPTION:
%   Invert a 1D warping function.
%
% INPUT:
%   gamma - A 1xN Matrix representing a 1D warping function
%
% OUTPUT:
%   gamma_i - The inverted warping function
%   Side effects: none
%
% ASSUMPTIONS AND LIMITATIONS:
%   None

function gamma_i = invert_gamma(gamma)

N = length(gamma);
x = (0:N-1) / (N-1);

gamma_i = interp1(gamma,x,x);

if isnan(gamma_i(N))
    gamma_i(N) = 1;
else
    gamma_i = gamma_i./gamma_i(N);
end

end
