function X = rand_von_mises(kappa, p, mu, N)

t = rand_marginal_t(kappa,p,N);
v = rand_uniform_hypersphere(N,p-1);

X = [repmat(sqrt(1 - t'.^2), p-1, 1).*v; t'];

pole = [zeros(p-1, 1); 1];
vecsum = mu + pole;

% Rodrigues formula to rotate from pole to mu
rot = 2 * ((vecsum*vecsum') / (vecsum'*vecsum)) - eye(p);

X = rot'*X;

end

function t = rand_marginal_t(kappa, p, N)

b = (p-1) / (2*kappa + sqrt(4*kappa^2 + (p-1)^2));
x0 = (1-b) / (1+b);
c = (kappa * x0) + ((p-1) * log(1 - x0^2));

% step 1 & step 2
nnow = N; 
last = 0;
t = zeros(N, 1);

while(true)
    ntrial = max(round(nnow*1.2), nnow+10) ;
    
    Z = betarnd((p-1)/2, (p-1)/2, ntrial, 1);
    U = rand(ntrial,1);
    W = (1 - (1+b)*Z) ./ (1 - (1-b)*Z);
        
    idx = (kappa * W) + ((p-1) * log(1 - x0*W)) - c >= log(U);
    
    first = last;
    last = min(N, first + sum(idx));
    
    buffer = W(idx);
    t((first+1):last) = buffer(1:(last-first));
    
    if last >= N
        break;
    else
        nnow = nnow-sum(idx);
    end
end

end
    
function v = rand_uniform_hypersphere(N, p)
% generate N uniformly distributed points on a (p-1)-dimensional hyper-sphere 

v = randn(p,N);
norms = zeros(N,1);
idx = ones(N,1,'logical');

% avoid vectors that are too small
% to prevent numerical discretisation
while any(idx)
    norms(idx) = vecnorm(v(:,idx));
    idx(idx) = norms(idx) < 1e-05;
   
    v(:,idx) = randn(p, sum(idx));
end

% normalise the vectors
v = v ./ norms';

end
    
%    References:
%
%    [1] Muller, M. E. "A Note on a Method for Generating Points Uniformly on N-Dimensional Spheres."
%    Comm. Assoc. Comput. Mach. 2, 19-20, Apr. 1959.
%
%    [2] https://mathworld.wolfram.com/SpherePointPicking.html
%
%    [3] https://dlwhittenbury.github.io/ds-2-sampling-and-visualising-the-von-mises-fisher-distribution-in-p-dimensions.html

