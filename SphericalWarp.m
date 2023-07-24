% SCRIPT NAME:
%   SphericalWarp
%
% DESCRIPTION:
%   Class to represent a diffeomorphism on the sphere.
%
% MATLAB VERSION:
%   R2022b
%
classdef SphericalWarp
    properties
        V
        J     
        T
    end
    properties (Access = private)        
        delta       
        P
        p_Vx 
        p_Vy 
        m_Vx 
        m_Vy
        p_Tx
        p_Ty 
        m_Tx 
        m_Ty               
       
        e1
        e2
        aabb        
    end
    methods
        function obj = SphericalWarp(grid, delta)
            obj.e1 = grid.e1;
            obj.e2 = grid.e2;
            obj.P = size(grid.V,1);

            % The warped vertices
            obj.V = grid.V;

            % Determinate of the Jacobian Matrix
            obj.J = ones(obj.P,1);       

            % The triangulation
            obj.T = grid.T;
     
            % Offset vertices along e1 and e2 (the tangent basis)
            p_Dx = sphere_exp_map(grid.V, grid.e1 * delta);
            p_Dy = sphere_exp_map(grid.V, grid.e2 * delta);
            m_Dx = sphere_exp_map(grid.V, grid.e1 * -delta);
            m_Dy = sphere_exp_map(grid.V, grid.e2 * -delta);
            
            obj.aabb = AABBtree(grid);
            
            % Interpolated values at delta coordinates
            [obj.p_Vx,obj.p_Tx] = obj.aabb.get_barycentric_data(p_Dx, false);
            [obj.m_Vx,obj.m_Tx] = obj.aabb.get_barycentric_data(m_Dx, false);
            [obj.p_Vy,obj.p_Ty] = obj.aabb.get_barycentric_data(p_Dy, false);
            [obj.m_Vy,obj.m_Ty] = obj.aabb.get_barycentric_data(m_Dy, false);

            % Small value for derivatives
            obj.delta = 2*delta;           
        end

        % Function to compose to diffeomorphisms
        function obj = compose_warp(obj, displacement)
            % Project displacement vector to the tangent plane in R3 then project to sphere 
            tangent_vec = (obj.e1 .* displacement(:,1) + obj.e2 .* displacement(:,2));            

            % Compose the new warp to the current overall warp
            [Vx,Tx] = obj.aabb.get_barycentric_data(obj.V, false);

            tangent_vec = obj.parallel_transport(tangent_vec(Tx,:)', obj.aabb.V(Tx,:)', repmat(obj.V,3,1)')';            
            tangent_vec = Vx(:) .* tangent_vec;
            tangent_vec = squeeze(sum(reshape(tangent_vec, obj.P, 3, 3), 2));
            
            % Project to sphere
            obj.V = normr(sphere_exp_map(obj.V, tangent_vec));      
            obj.J = obj.get_jacobian();
        end        

	function obj = invert_warp(obj)
            tangent_vec = sphere_log_map(obj.V,obj.aabb.V);

            grid.V = obj.V;
            grid.T = obj.T;

            tmp_aabb = AABBtree(grid);

            % Compose the new warp to the current overall warp
            [Vx,Tx] = tmp_aabb.get_barycentric_data(obj.aabb.V, false);

            tangent_vec = obj.parallel_transport(tangent_vec(Tx,:)', obj.V(Tx,:)', repmat(obj.aabb.V,3,1)')';
            tangent_vec = Vx(:) .* tangent_vec;
            tangent_vec = squeeze(sum(reshape(tangent_vec, obj.P, 3, 3), 2));

            % Project to sphere
            obj.V = normr(sphere_exp_map(obj.V, tangent_vec));
            obj.J = obj.get_jacobian();
        end

        function plot(obj, fig_title)            
            title(fig_title)
            trisurf(obj.T, obj.V(:,1), obj.V(:,2), obj.V(:,3), obj.J)
            axis off
            axis equal
        end
    end

    methods (Access = private)
        % Function to compute the Jacobian of a diffeomorphism
        function J = get_jacobian(obj)
            % Interpolate warp in positive theta direction
            wx = dot(obj.p_Vx, reshape(obj.V(obj.p_Tx,1),obj.P,3),2);
            wy = dot(obj.p_Vx, reshape(obj.V(obj.p_Tx,2),obj.P,3),2);
            wz = dot(obj.p_Vx, reshape(obj.V(obj.p_Tx,3),obj.P,3),2);            

            [p_Dtt,p_Dpt] = cart_to_sphere(normr([wx,wy,wz]));

            % Interpolate warp in negative theta direction
            wx = dot(obj.m_Vx, reshape(obj.V(obj.m_Tx,1),obj.P,3),2);
            wy = dot(obj.m_Vx, reshape(obj.V(obj.m_Tx,2),obj.P,3),2);
            wz = dot(obj.m_Vx, reshape(obj.V(obj.m_Tx,3),obj.P,3),2);            

            [m_Dtt,m_Dpt] = cart_to_sphere(normr([wx,wy,wz]));

            % Interpolate warp in positive phi direction
            wx = dot(obj.p_Vy, reshape(obj.V(obj.p_Ty,1),obj.P,3),2);
            wy = dot(obj.p_Vy, reshape(obj.V(obj.p_Ty,2),obj.P,3),2);
            wz = dot(obj.p_Vy, reshape(obj.V(obj.p_Ty,3),obj.P,3),2);            

            [p_Dtp,p_Dpp] = cart_to_sphere(normr([wx,wy,wz]));

            % Interpolate warp in negative phi direction
            wx = dot(obj.m_Vy, reshape(obj.V(obj.m_Ty,1),obj.P,3),2);
            wy = dot(obj.m_Vy, reshape(obj.V(obj.m_Ty,2),obj.P,3),2);
            wz = dot(obj.m_Vy, reshape(obj.V(obj.m_Ty,3),obj.P,3),2);            

            [m_Dtp,m_Dpp] = cart_to_sphere(normr([wx,wy,wz]));  

            % Deal with the boundry cases of theta,phi by finding the
            % minimum absolute value of the difference +- j*pi
            offset = ((-2:2).*pi);

            % The Jacobian matrix [Dtt Dtp; Dpt Dpp] using finite difference, t=theta, p=phi
            Dtt = min(p_Dtt - m_Dtt + offset,[],2,'ComparisonMethod','abs') ./ (obj.delta);
            Dtp = min(p_Dtp - m_Dtp + offset,[],2,'ComparisonMethod','abs') ./ (obj.delta);
            Dpt = min(p_Dpt - m_Dpt + offset,[],2,'ComparisonMethod','abs') ./ (obj.delta);
            Dpp = min(p_Dpp - m_Dpp + offset,[],2,'ComparisonMethod','abs') ./ (obj.delta);

            % The determinant of the Jacobian
            J = ((Dtt.*Dpp) - (Dtp.*Dpt));

            % The Jacobian operation is a linear mapping between the tangent
            % spaces T_(t,p) to T_(wt,wp). So to get the Jacobian matrix with
            % respect to the sphere at (t,p), we need to apply a change of
            % basis to the linear Jacobian matrix. First we apply a mapping
            % to transform from T_(wt,wp) into e1,e2, and then get the
            % Jacobian. Finally, we transform back to elements of T_(t,p).
            % In matrix form [1 0;1 sin(wt)] J_(t,p) [1 0; 0 1/sin(t)]
            [theta,~] = cart_to_sphere(obj.V);
            J = abs((J .* (sin(theta))));
        end

        function new_vec = parallel_transport(~, tan_vec, orig_pt, new_pt)
            tol = 1e-3;

            % Ensure that tan_vec is tangent to the orig_pts
            zero_vec = max(abs(dot(tan_vec, orig_pt, 1) ./ (vecnorm(tan_vec) + eps)));
            
            if(zero_vec > tol)
                error('tan_vec must be perpendicular to orig_pt');                
            end

            % rotation axis between orig_pt and new_pt
            rotation = cross(orig_pt, new_pt, 1);
            rotation = rotation ./ (vecnorm(rotation) + eps);

            % angle between rotation axis and orig_pt
            old_v = cross(rotation, orig_pt, 1);
            old_v = old_v ./ (vecnorm(old_v) + eps);

            % angle between rotation axis and new_pt
            new_v = cross(rotation, new_pt, 1);
            new_v = new_v ./ (vecnorm(new_v) + eps);

            % project tan_vec onto new points
            proj_on_newv = dot(tan_vec, old_v, 1);
            proj_on_neww = dot(tan_vec, rotation, 1);

            % the transported vector
            new_vec = repmat(proj_on_newv, 3, 1) .* new_v + repmat(proj_on_neww, 3, 1) .* rotation;

            % clean up close points
            idx = (vecnorm(orig_pt - new_pt) < 1e-4);
            new_vec(:, idx) = tan_vec(:, idx);

            % clean up antipodal points 
            idx = (vecnorm(orig_pt + new_pt) < 1e-4);
            new_vec(:, idx) = tan_vec(:, idx);
        end
    end
end
