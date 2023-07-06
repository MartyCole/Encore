classdef SphericalDerivative    
    properties (Access = private)        
        delta
        delta2
        P
        p_Vx 
        p_Vy 
        m_Vx 
        m_Vy
        p_Tx
        p_Ty 
        m_Tx 
        m_Ty 
        dVx
        dTx
        dVy
        dTy
        Vg
        Tg
        direction
        interpolator
        grid
    end
    methods
        function obj = SphericalDerivative(grid, delta, direction)          
            % Small value for derivatives
            obj.delta  = delta;
            obj.delta2 = delta*2;

            % Grid theta, phi, and size P
            [theta,phi] = obj.cart_to_sphere(grid.V); 
            obj.grid = grid;
            obj.P = length(grid.V);

            obj.direction = direction;
            if strcmp(direction, 'polar')
                % Offset vertices along theta and phi
                p_Dx = sphere_to_cart(theta+delta,phi).';
                p_Dy = sphere_to_cart(theta,phi+delta).';
                m_Dx = sphere_to_cart(theta-delta,phi).';
                m_Dy = sphere_to_cart(theta,phi-delta).';
            elseif strcmp(direction, 'tangent')
                % Offset vertices along e1 and e2 (the tangent basis)
                p_Dx = sphere_exp_map(grid.V, grid.e1 * delta);
                p_Dy = sphere_exp_map(grid.V, grid.e2 * delta);
                m_Dx = sphere_exp_map(grid.V, grid.e1 * -delta);
                m_Dy = sphere_exp_map(grid.V, grid.e2 * -delta);
            else
                error('Direction must be "polar", or "tangent"')                
            end

            obj.interpolator = SphericalInterpolator(grid);
            
            % Interpolated values at delta coordinates
            [obj.p_Vx,obj.p_Tx] = obj.interpolator.get_barycentric_data(p_Dx);
            [obj.m_Vx,obj.m_Tx] = obj.interpolator.get_barycentric_data(m_Dx);
            [obj.p_Vy,obj.p_Ty] = obj.interpolator.get_barycentric_data(p_Dy);
            [obj.m_Vy,obj.m_Ty] = obj.interpolator.get_barycentric_data(m_Dy);

            obj.dVx = obj.p_Vx.';
            obj.dTx = int64(obj.p_Tx-1).';
            obj.dVy = obj.p_Vy.';
            obj.dTy = int64(obj.p_Ty-1).';

            [obj.Vg,obj.Tg] = obj.interpolator.get_barycentric_data(grid.V);
            obj.Vg = obj.Vg.';
            obj.Tg = int64(obj.Tg-1).';
        end

        function [dQxe1,dQxe2,dQye1,dQye2] = get_derivative(obj, Q)
            dQxe1 = (bary_interp_2D_mex(obj.dVx,obj.Vg,obj.dTx,obj.Tg,Q) - Q) / obj.delta;
            dQxe2 = (bary_interp_2D_mex(obj.dVy,obj.Vg,obj.dTy,obj.Tg,Q) - Q) / obj.delta;
            dQye1 = (bary_interp_2D_mex(obj.Vg,obj.dVx,obj.Tg,obj.dTx,Q) - Q) / obj.delta;
            dQye2 = (bary_interp_2D_mex(obj.Vg,obj.dVy,obj.Tg,obj.dTy,Q) - Q) / obj.delta;
        end     

        % Function to compose to diffeomorphisms
        function new_warp = compose_warp(obj, displacement, warp)
            % Project displacement vector to the tangent plane in R3 then project to sphere 
            tangent_vec = (obj.grid.e1 .* displacement(:,1) + obj.grid.e2 .* displacement(:,2));
            new_points = normr(sphere_exp_map(obj.grid.V, tangent_vec)); 

            % Compose the new warp to the current overall warp
            [Vx,Tx] = obj.interpolator.get_barycentric_data(new_points);
            
            v1 = warp.V(Tx(:,1),:);
            v2 = warp.V(Tx(:,2),:);
            v3 = warp.V(Tx(:,3),:);

            % Project the interpolated displacement vector to the sphere
            new_warp = warp;
            new_warp.V = normr(Vx(:,1).*v1 + Vx(:,2).*v2 + Vx(:,3).*v3);         
        end

        % Function to compute the Jacobian of a diffeomorphism
        function J = get_jacobian(obj, warp, map_to_T)
            % Interpolate warp in positive theta direction
            wx = dot(obj.p_Vx, reshape(warp.V(obj.p_Tx,1),obj.P,3),2);
            wy = dot(obj.p_Vx, reshape(warp.V(obj.p_Tx,2),obj.P,3),2);
            wz = dot(obj.p_Vx, reshape(warp.V(obj.p_Tx,3),obj.P,3),2);            

            [p_Dtt,p_Dpt] = obj.cart_to_sphere(normr([wx,wy,wz]));

            % Interpolate warp in negative theta direction
            wx = dot(obj.m_Vx, reshape(warp.V(obj.m_Tx,1),obj.P,3),2);
            wy = dot(obj.m_Vx, reshape(warp.V(obj.m_Tx,2),obj.P,3),2);
            wz = dot(obj.m_Vx, reshape(warp.V(obj.m_Tx,3),obj.P,3),2);            

            [m_Dtt,m_Dpt] = obj.cart_to_sphere(normr([wx,wy,wz]));

            % Interpolate warp in positive phi direction
            wx = dot(obj.p_Vy, reshape(warp.V(obj.p_Ty,1),obj.P,3),2);
            wy = dot(obj.p_Vy, reshape(warp.V(obj.p_Ty,2),obj.P,3),2);
            wz = dot(obj.p_Vy, reshape(warp.V(obj.p_Ty,3),obj.P,3),2);            

            [p_Dtp,p_Dpp] = obj.cart_to_sphere(normr([wx,wy,wz]));

            % Interpolate warp in negative phi direction
            wx = dot(obj.m_Vy, reshape(warp.V(obj.m_Ty,1),obj.P,3),2);
            wy = dot(obj.m_Vy, reshape(warp.V(obj.m_Ty,2),obj.P,3),2);
            wz = dot(obj.m_Vy, reshape(warp.V(obj.m_Ty,3),obj.P,3),2);            

            [m_Dtp,m_Dpp] = obj.cart_to_sphere(normr([wx,wy,wz]));  

            % Deal with the boundry cases of theta,phi by finding the
            % minimum absolute value of the difference +- j*pi
            offset = ((-2:2).*pi);

            % The Jacobian matrix [Dtt Dtp; Dpt Dpp] using finite difference, t=theta, p=phi
            Dtt = min(m_Dtt - p_Dtt + offset,[],2,'ComparisonMethod','abs') ./ (obj.delta2);
            Dtp = min(m_Dtp - p_Dtp + offset,[],2,'ComparisonMethod','abs') ./ (obj.delta2);
            Dpt = min(m_Dpt - p_Dpt + offset,[],2,'ComparisonMethod','abs') ./ (obj.delta2);
            Dpp = min(m_Dpp - p_Dpp + offset,[],2,'ComparisonMethod','abs') ./ (obj.delta2);            

            % The determinant of the Jacobian
            J = ((Dtt.*Dpp) - (Dtp.*Dpt));

            % The Jacobian operation is a linear mapping between the tangent 
            % spaces T_(t,p) to T_(wt,wp). So to get the Jacobian matrix with 
            % respect to the sphere at (t,p), we need to apply a change of 
            % basis to the linear Jacobian matrix. First we apply a mapping
            % to transform from T_(wt,wp) into e1,e2, and then get the
            % Jacobian. Finally, we transform back to elements of T_(t,p).
            % In matrix form [1 0;1 sin(wt)] J_(t,p) [1 0; 0 1/sin(t)]            
            if map_to_T 
                [theta,~] = obj.cart_to_sphere(warp.V);

                if strcmp(obj.direction, 'polar')
                    J = (J .* (sin(theta)./sin(obj.grid.theta)));
                else
                    J = (J .* (sin(theta)));
                end
            end            
               
            % The final value
            J = abs(J);            
        end
    end
    methods (Access = private)
        function [theta, phi] = cart_to_sphere(~, x)
            theta = acos(x(:,3));
            phi = atan2(x(:,2),x(:,1));

            phi(phi < 0) = phi(phi < 0) + 2*pi;  
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