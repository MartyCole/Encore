function [ddPlmdx, dPlmdx, Plm] = legendre_2nd_derivative(l,x)

    %% Initialization and checks
    assert(isnumeric(l) && isscalar(l) && isfinite(l) && isreal(l) && round(l)==l && l>=0,...
          [mfilename ':invalid_l'],...
          'Degree L must be a positive scalar integer.');
      
    assert(~isempty(x) && isnumeric(x) && all(isreal(x(:))) && ismatrix(x),...
           [mfilename ':invalid_x'],...
           'X must be real, numeric, nonempty, and be a vector or matrix.');
    
    % Make sure that if x is a 
    % vector, it is a row vector
    if size(x,2) == 1
        x = x';
    end

    m = (0:l)';

    % Normalisation constant for Spherical Harmonics
    Nlm = sqrt(((2*l+1) / (4*pi)) * (factorial(l-m) ./ factorial(l+m)));
    Nlminus = sqrt((l+m).*(l-m+1));
    Nlplus  = sqrt((l-m).*(l+m+1));
        
    % Calculate Plm so that we can use a
    % recursion relation to calculate dPlmdx  
    Plm = Nlm .* legendre(l,cos(x));
                     
    %% Computation 
    if numel(Plm)==1
        dPlmdx = zeros(size(x)); 
        return; 
    end    
  
    % Special case of vector Plm
    if size(x,1) == 1
        Plm = permute(Plm, [1 3 2]); 
    end     

    % Calculate Pln for n = [-1,0,...,m-1]
    Plminus = [ -Plm(2,:,:) / (l*(l+1)) 
           Plm(1:end-1,:,:)];
       
    % Calculate Pln for n = [1,...,m+1]
    Plplus = [ Plm(2:end,:,:); zeros(1,size(Plm,2),size(Plm,3))];     

    % Calculate the derivative
    m = (0:l)';
    dPlmdx = -0.5 * ((Plminus .* Nlminus) - (Plplus .* Nlplus));
    
    % Calculate dPlndx for n = [-1,0,...,m-1]
    dPlminus = [ dPlmdx(2,:,:) / (l*(l+1)) 
           -dPlmdx(1:end-1,:,:)];
       
    % Calculate dPlndx for n = [1,...,m+1]
    dPlplus = [ dPlmdx(2:end,:,:); zeros(1,size(Plm,2),size(Plm,3))];     
    
    % Calculate second derivative
    ddPlmdx = 0.5 * ((dPlminus .* Nlminus) - (dPlplus .* Nlplus));    
end