classdef SphericalGrid
    % SphericalGrid
    %
    % This class represents a spherical discretization used for registration
    % and kernel evaluation on the sphere. Given a triangular mesh approximating
    % a sphere, the class computes:
    %
    %   - normalized vertex positions on the unit sphere
    %   - spherical coordinates (theta, phi) for each vertex
    %   - Voronoi areas for numerical integration
    %   - a complete tangent‑space basis system derived from spherical harmonics
    %   - divergence of each basis field
    %   - an orthonormal Cartesian tangent frame (e1, e2) at every vertex
    %
    % The tangent basis consists of:
    %   - gradient fields of real spherical harmonics
    %   - curl fields of the same harmonics
    %   - area‑normalized basis vectors suitable for variational optimization
    %
    % The class also provides private utilities for:
    %   - computing spherical harmonics and their derivatives
    %   - constructing gradient and curl fields
    %   - normalizing basis functions using area‑weighted L2 norms
    %
    % This grid serves as the geometric foundation for the Encore registration 
    % algorithms, enabling differential operators, kernel evaluation, and
    % tangent‑space optimization on the sphere.
    %
    properties
        V % Vertex positions (cartesian)
        T % Triangle (face) indices
        A % Voronoi area around each vertex

        theta % Polar angle (colatitude)
        phi % Azimuth angle   

        num_basis % Number of basis functions (2(L+1)^2 - 2)
        basis % Tangent basis field vectors
        divergence % Divergence of each basis field

        e1 % First tangent vector (theta direction)
        e2 % Second tangent vector (phi direction)

        use_GPU logical = true % convert matrices to GPU when enabled
    end 
    methods
        function obj = SphericalGrid(mesh, L, gpu)
            % Construct a spherical grid from a triangular mesh.
            %
            % Inputs:
            %   mesh : struct with fields V (vertices) and T (triangles)
            %   L    : maximum spherical harmonic degree used for basis construction

            if nargin < 3
                gpu = true;
            end
            
            obj.use_GPU = gpu;

            % normalise mesh vertices to lie on the unit sphere
            obj.V = normr(mesh.V);
            obj.T = mesh.T;
            
            % convert from Cartesian to spherical coordinates
            [obj.theta, obj.phi] = cart_to_sphere(obj.V);

            % compute Voronoi area for each vertex (used for integration)
            obj.A = calc_voronoi_area(obj.V, obj.T);
            
            % build tangent-space basis fields (gradients + curls)
            [basis, divergence] = obj.build_tangent_basis(L);

            % move basis fields to GPU for faster computation
            obj.basis = obj.maybe_GPU(basis);
            obj.divergence = obj.maybe_GPU(divergence);

            % store the number of basis functions
            obj.num_basis = size(obj.basis,2);

            % construct orthonormal tangent frame (e1, e2) in Cartesian 
            % coordinates, where e1 = d/dtheta, e2 = d/dphi
            obj.e1 = obj.maybe_GPU([cos(obj.theta).*cos(obj.phi), ...
                               cos(obj.theta).*sin(obj.phi), ...
                              -sin(obj.theta)]);
    
            obj.e2 = obj.maybe_GPU([-sin(obj.phi), ...
                                cos(obj.phi), ...
                                zeros(size(obj.theta,1),1)]);
        end
    end

    methods (Access = private) 
        function [basis, div] = build_tangent_basis(obj, L)
            % Build tangent basis fields from spherical harmonics up to degree L.
            %
            % For each harmonic Y_lm:
            %   grad(Y) = [dY/dtheta, dY/dphi]
            %   curl(Y) = [dY/dphi, -dY/dtheta]
            %
            % Divergence:
            %   div(grad(Y)) = Laplacian(Y)
            %   div(curl(Y)) = 0
            %
            % Steps:
            %   - compute derivatives and Laplacian of harmonics for each l = 1..L
            %   - assemble gradient basis fields
            %   - normalize basis fields using area‑weighted L2 norm
            %   - append curl basis fields
            %   - compute divergence for each basis
            %
            % Outputs:
            %   basis : NxBx2 array of tangent basis vectors
            %   div   : NxB divergence values

            K = (L+1)^2 - 1;          
            N = size(obj.theta, 1);        

            % allocate basis arrays
            basis = zeros(N, 2*K, 2);
            div   = zeros(N, 2*K);

            idx = 1;

            % loop over harmonic degree l = 1..L
            for l = 1:L
                % number of modes, m, for degree l
                len = 2*(l+1) - 1;

                % compute derivatives and Laplacian of harmonics
                [dY_theta, dY_phi, lap_Y] = obj.compute_harmonics(l);
                
                % store gradient components
                basis(:,idx:(idx+len-1),1) = dY_theta;
                basis(:,idx:(idx+len-1),2) = dY_phi;

                % store divergence (Laplacian)
                div(:,idx:(idx+len-1)) = lap_Y;

                idx = idx + len;
            end

            % normalize basis fields using area-weighted L2 norm
            norm_constant = sqrt(sum((basis(:,:,1).^2 + basis(:,:,2).^2) .* obj.A, 1));            
                 
            basis(:,:,1) = basis(:,:,1) ./ norm_constant;
            basis(:,:,2) = basis(:,:,2) ./ norm_constant;     
            div = div ./ norm_constant;        
            
            % build the curl basis fields: curl(Y) = [dY/dphi, -dY/dtheta]
            basis(:,(K+1):end,1) = basis(:,1:K,2);
            basis(:,(K+1):end,2) = -basis(:,1:K,1);

            % divergence of curl fields is zero
            div(:,(K+1):end) = 0;
        end

        function [dY_theta, dY_phi, lap_Y] = compute_harmonics(obj, L)    
            % Compute real spherical harmonics of degree L and their derivatives.
            %
            % Outputs:
            %   dY_theta : derivative of Y with respect to theta
            %   dY_phi   : derivative of Y with respect to phi
            %   lap_Y    : Laplacian of Y (−L(L+1)Y)
            %
            % Internally computes:
            %   - normalized associated Legendre polynomials P_lm
            %   - stable derivatives using recurrence relations
            %   - real harmonic indexing: Re(Y_0), Re(Y_1), Im(Y_1), Re(Y_2), ...

            % compute normalised associated Legendre polynomials and derivatives
            [dPlm_dx, Plm] = obj.legendre_derivatives(L, cos(obj.theta));

            P = 2*(L+1) - 1;            
            m = (0:L);
            N = size(obj.theta,1);
    
            % precompute trig terms
            sin_phi = (sin(obj.phi * m));
            cos_phi = (cos(obj.phi * m));
            sin_theta = (max(sin(obj.theta), 1e-8)); % avoid division by xero
            
            % allocate outputs
            Y = zeros(N, P);   
            dY_theta = zeros(N, P);            
            dY_phi = zeros(N, P);
            
            % indexing for real/imag parts: Re(Y_0), Re(Y_1), Im(Y_1), Re(Y_2), Im(Y_2), ...
            real_idx = [1, 2:2:P];
            im_idx = 3:2:P;

            % real parts: cos(m*phi)
            Y(:,real_idx) = Plm .* cos_phi; 
            dY_theta(:,real_idx) = dPlm_dx .* cos_phi; 
            dY_phi(:,real_idx) = Plm .* -(m .* sin_phi) ./ sin_theta; 
            
            % imaginary parts: sin(m*phi), m > 0 
            Y(:,im_idx) = Plm(:,2:end) .* sin_phi(:,2:end); 
            dY_theta(:,im_idx) = dPlm_dx(:,2:end) .* sin_phi(:,2:end); 
            dY_phi(:,im_idx) = Plm(:,2:end) .* (m(2:end) .* cos_phi(:,2:end)) ./ sin_theta;   
            
            lap_Y = -L*(L + 1) * Y;            
        end   

        function [dPlm_dx, Plm] = legendre_derivatives(~, L, X)   
            % Compute normalized associated Legendre polynomials P_lm and their
            % derivatives with respect to x = cos(theta).
            %
            % Uses stable recurrence relations:
            %   - P_{l,m−1} and P_{l,m+1}
            %   - derivative formula involving Nlm_minus and Nlm_plus
            %
            % Outputs:
            %   dPlm_dx : derivative of P_lm with respect to x
            %   Plm     : normalized associated Legendre polynomials

            N = size(X,1);
            m = (0:L);
        
            % normalisation constant for spherical harmonics
            Nl = sqrt(((2*L+1) / (4*pi)) * (factorial(L-m) ./ factorial(L+m))).';
            Nlm_minus = sqrt((L+m).*(L-m+1));
            Nlm_plus  = sqrt((L-m).*(L+m+1));
                
            % compute normalized associated Legendre polynomials
            Plm = (Nl .* legendre(L,X)).';                             
          
            % calculate P_{l,m-1} using recurrence
            Plm_minus = [ -Plm(:,2)/(L*(L+1)), Plm(:,1:end-1) ];            
               
            % calculate P_{l,m+1}
            Plm_plus = [ Plm(:,2:end), zeros(N,1) ]; 
        
            % calculate the derivative using stable recurrence relation
            dPlm_dx = -0.5 * (Plm_minus .* Nlm_minus - Plm_plus .* Nlm_plus);
        end

        function X = maybe_GPU(obj, X)
            % Convert matrix to GPU if enabled.            
            %
            % Output:
            %   X : possibly GPU‑resident matrix

            if obj.use_GPU == true
                X = gpuArray(X);
            end
        end
    end
end