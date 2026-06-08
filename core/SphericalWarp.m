classdef SphericalWarp < matlab.mixin.Copyable
    % SphericalWarp
    %
    % This class represents a diffeomorphic warp (smooth, invertible mapping)
    % on the sphere. A warp is defined by a tangent‑plane displacement field
    % at each vertex (a stationary velocity field, SVF), and is updated using:
    %   - the spherical exponential map (scale and squaring)
    %   - parallel transport of tangent vectors
    %   - barycentric interpolation via an AABB tree
    %
    % The class supports:
    %   - composing incremental warp updates
    %   - computing the inverse warp
    %   - applying rigid rotations
    %   - visualizing the warp field and its Jacobian determinant
    %
    % Used by Encore to apply and update diffeomorphic transformations
    % during spherical registration.
    %

    properties
        v % stationary velocity field
        V % warped vertex positions
        T % triangulation
    end
    properties (Dependent)
        J % Jacobian determinant (area distortion)
    end
    properties (Access = private)    
        P % number of vertices       
       
        V0 % reference (unwarped) vertices
        L % Laplacian for smoothing the SVF
        e1 % tangent frame vector
        e2 % tangent frame vector
        grid_A % original Voronoi areas (for Jacobian computation)
    end
    properties (Access = private, NonCopyable)
        aabb % AABB tree for barycentric interpolation    
    end
    methods
        function obj = SphericalWarp(grid)
            % Construct a SphericalWarp object from a spherical grid.
            %
            % Initializes:
            %   - tangent frame (e1, e2)
            %   - vertex positions and triangulation
            %   - original Voronoi areas for Jacobian computation
            %   - Laplacian for smoothing the SVF
            %   - AABB tree for barycentric interpolation
            %
            % Input:
            %   grid : SphericalGrid object defining the base mesh

            obj.e1 = grid.e1;
            obj.e2 = grid.e2;
            obj.P = size(grid.V,1);

            % the warped vertices
            obj.V = grid.V;                        

            % the triangulation
            obj.T = grid.T;            

            obj.grid_A = calc_voronoi_area(obj.V, obj.T);            
            obj.aabb = AABBtree(grid);

            % SVF displacement c and reference vertices V0
            obj.v = zeros(obj.P, 2);
            obj.V0 = obj.V;

            % Laplacian for smoothing
            obj.L = obj.cot_matrix();
        end
        
        function obj = compose(obj, displacement)
            % Compose an incremental displacement field with the current warp.
            %
            % Steps:
            %   - add displacement to the stationary velocity field
            %   - smooth the field to keep things well‑behaved
            %   - apply the exponential map to get updated vertex positions
            %   - reject updates that would create flipped triangles
            %
            % Input:
            %   displacement : Px2 tangent displacement field

            % composition is just addition in SVFs
            obj.v = obj.v + displacement;

            % apply some "viscous-ness" to the field to keep it smooth
            obj.v = obj.v - 0.05 * (obj.L * double(obj.v));

            % apply the exponention map to get vertex positions
            test_V = obj.exp_svf();

            % don't apply warp if it results in folded triangles
            p1 = test_V(obj.T(:,1),:)';
            p2 = test_V(obj.T(:,2),:)';
            p3 = test_V(obj.T(:,3),:)';

            normals = dot(cross(p2 - p1, p3 - p1), p1+p2+p3);
            signs = sign(normals);

            if any(signs <= 0)
                fprintf("WARNING: Skipping update to avoid folded triangles.\n")
            else
                obj.V = test_V;
            end
        end
        
	    function obj = invert(obj)
            % Compute the inverse of the current warp.
            %
            % Steps:
            %   - negate the stationary velocity field gives the inverse flow
            %   - apply the exponential map to produce the inverse warp         
            %
            % Produces a warp such that: inverse_warp o warp ~ identity

            obj.v = -obj.v;
            obj.V = obj.exp_svf();      
        end

        function obj = resample(obj, new_grid)
            % tangent velocity in euclidean space
            tangent_vec = obj.e1 .* obj.v(:,1) + obj.e2 .* obj.v(:,2);  

            % interpolate 3D velocity onto new grid vertices
            [Vx, Tx] = obj.aabb.get_barycentric_data(new_grid.V, false);

            v_interp = Vx(:) .* tangent_vec(Tx,:);
            v_interp = squeeze(sum(reshape(v_interp, size(new_grid.V,1), 3, 3), 2)); 

            % project interpolated velocity into new tangent frames
            obj.v = [dot(v_interp, new_grid.e1, 2), dot(v_interp, new_grid.e2, 2)];

            % update reference geometry
            obj.V0 = new_grid.V;
            obj.e1 = new_grid.e1;
            obj.e2 = new_grid.e2;
            obj.T  = new_grid.T;
            obj.P  = size(new_grid.V,1);
            obj.aabb = AABBtree(new_grid);
           
            % recompute warped vertices from SVF
            obj.V = obj.exp_svf(obj.v, obj.V0);
        end

        function obj = rotate(obj, R)
            % Apply a rigid rotation to the warp.
            %
            % Steps:
            %   - converts the rotation matrix to an axis‑angle vector
            %   - compute the induced tangent velocity field on the sphere
            %   - update the warp via the exponential map
            %
            % Input:
            %   R : 3x3 rotation matrix (SO(3))           
     
            theta = acos(max(-1, min(1, (trace(R) - 1) / 2)));
        
            if abs(theta) < 1e-12
                omega = [0 0 0];                
            else 
                wx = (R(3,2) - R(2,3)) / (2*sin(theta));
                wy = (R(1,3) - R(3,1)) / (2*sin(theta));
                wz = (R(2,1) - R(1,2)) / (2*sin(theta));
            
                axis = [wx wy wz];
            
                % rotation vector
                omega = axis * theta;       
            end

            % 3D velocity field on the sphere
            v3 = cross(repmat(omega, size(obj.V0,1), 1), obj.V0, 2);
            
            % project to tangent coordinates
            v1 = dot(v3, obj.e1, 2);
            v2 = dot(v3, obj.e2, 2);
            
            obj.v = [v1, v2];   

            % map onto vertices
            obj.V = obj.exp_svf();
        end

        function plot(obj, fig_title) 
            % Visualize the warped surface colored by its Jacobian determinant.
            %
            % Inputs:
            %   fig_title : title for the figure
            %
            % Displays:
            %   - trisurf of warped mesh
            %   - color map representing local area distortion (J)

            trisurf(obj.T, obj.V(:,1), obj.V(:,2), obj.V(:,3), obj.J)
            title(fig_title)
            axis off
            axis equal
        end

        function plot_field(obj, fig_title)
            % Visualize the displacement field of the warp.
            %
            % Displays:
            %   - quiver3 plot of tangent vectors (log map of warp)
            %   - underlying warped mesh for context
            %
            % Input:
            %   fig_title : title for the figure

            tangent_vec = obj.e1 .* obj.v(:,1) + obj.e2 .* obj.v(:,2);
            mag = sqrt(sum(tangent_vec.^2, 2));

            % Normalize colors to [0,1]
            mag_norm = (mag - min(mag)) / (max(mag) - min(mag) + eps);

            quiver3(obj.V0(:,1),obj.V0(:,2),obj.V0(:,3), ...
                tangent_vec(:,1), tangent_vec(:,2), tangent_vec(:,3),2,'k','LineWidth',2);

            hold on; axis equal; colormap jet; colorbar;

            trisurf(obj.T, obj.V0(:,1)*0.99, obj.V0(:,2)*0.99, obj.V0(:,3)*0.99, mag_norm, ...
                'EdgeAlpha', 0)

            shading interp;

            hold off

            title(fig_title)            
        end    

        function J = get.J(obj)    
            % Compute the Jacobian determinant of the warp.
            %
            % J = (warped Voronoi area) / (original Voronoi area)
            %
            % Represents local area distortion induced by the diffeomorphism.

            A1 = calc_voronoi_area(obj.V, obj.T);

            % Mean area ratio = Jacobian determinant
            J = (A1 ./ obj.grid_A);
        end
    end

    methods (Access = private)
        function V = exp_svf(obj)
            % Apply the exponential map to the stationary velocity field.
            %
            % Uses scaling‑and‑squaring:
            %   - scale the velocity field
            %   - repeatedly interpolate and parallel‑transport tangent vectors
            %   - apply the spherical exponential map at each squaring step
            %
            % Returns:
            %   V : updated vertex positions on the sphere

            % number of squaring steps
            N = 6;  
          
            % scaling step
            v_scaled = obj.v / (2^N);
          
            % perform the initial exponential step
            tangent_vec = obj.e1 .* v_scaled(:,1) + obj.e2 .* v_scaled(:,2);
            V = normr(sphere_exp_map(obj.V0, tangent_vec));

            % repeatedly square the deformation
            for i = 1:N              
                [Vx, Tx] = obj.aabb.get_barycentric_data(V, false);

                % compute the current tangent displacement wrt grid
                tan_vec_full = sphere_log_map(obj.V0, V);            
                
                % transport the displacement to the deformed sphere
                tan_vec_interp = obj.parallel_transport(tan_vec_full(Tx,:), obj.V0(Tx,:), repmat(V,3,1));

                % interpolate the SVF at the points on the deformed sphere
                tan_vec_interp = Vx(:) .* tan_vec_interp;
                tan_vec_interp = squeeze(sum(reshape(tan_vec_interp, obj.P, 3, 3), 2));

                % apply the exponential map (squaring step)
                V = normr(sphere_exp_map(V, tan_vec_interp));
            end
        end

        function new_vec = parallel_transport(~, tan_vec, orig_pt, new_pt)
            % Parallel transport a tangent vector from one point on the sphere to another.
            %
            % Uses the closed‑form formula for parallel transport on S^2:
            %   v_new = v - ((v·new_pt)/(1 + orig_pt·new_pt)) * (orig_pt + new_pt)
            %
            % Handles near‑antipodal cases by falling back to the original vector.
            %
            % Inputs:
            %   tan_vec : tangent vectors to transport
            %   orig_pt : original base points
            %   new_pt  : destination base points
            %
            % Output:
            %   new_vec : transported tangent vectors

            % dot products
            pq = dot(orig_pt, new_pt, 2);
            vq = dot(tan_vec, new_pt, 2);

            denom = 1 + pq;

            % compute transported vector
            correction = (vq ./ denom) .* (orig_pt + new_pt);
            new_vec = tan_vec - correction;

            % handle antipodal or near-antipodal cases
            idx_bad = denom < 1e-8;
            new_vec(idx_bad,:) = tan_vec(idx_bad,:);
        end

        function L = cot_matrix(obj)
            % Cotangent Laplacian for a triangular mesh

            i1 = obj.T(:,1); 
            i2 = obj.T(:,2); 
            i3 = obj.T(:,3);

            % edge vectors
            v1 = obj.V0(i2,:) - obj.V0(i3,:);
            v2 = obj.V0(i3,:) - obj.V0(i1,:);
            v3 = obj.V0(i1,:) - obj.V0(i2,:);

            % triangle areas
            n  = cross(v1, -v3, 2);
            dblA = sqrt(sum(n.^2, 2));

            % cotangents of angles
            cot12 = dot(v1, -v3, 2) ./ dblA; 
            cot23 = dot(v2, -v1, 2) ./ dblA;
            cot31 = dot(v3, -v2, 2) ./ dblA;  

            % build symmetric weights for edges
            I = [i2; i3; i3; i1; i1; i2];
            K = [i3; i2; i1; i3; i2; i1];
            S = [cot12; cot12; cot23; cot23; cot31; cot31] * 0.5;

            W = sparse(I, K, S, obj.P, obj.P);

            % Laplacian: L = D - W
            d = sum(W,2);
            L = spdiags(d, 0, obj.P, obj.P) - W;
        end
    end

    methods(Access = protected)
        function cpObj = copyElement(obj)
            % Deep copy method for SphericalWarp.
            %
            % Ensures that the AABB tree is reconstructed rather than shallow‑copied.
            %
            % Called automatically when using the Copyable mixin.

            cpObj = copyElement@matlab.mixin.Copyable(obj);

            grid.V = obj.aabb.V;
            grid.T = obj.aabb.T;

            cpObj.aabb = AABBtree(grid);
        end
    end

    % Hide unneeded documentation
    methods (Hidden)
        function varargout = findobj(O,varargin)
            varargout = findobj@handle(O,varargin{:});
        end
        function varargout = findprop(O,varargin)
            varargout = findprop@handle(O,varargin{:});
        end
        function varargout = addlistener(O,varargin)
            varargout = addlistener@handle(O,varargin{:});
        end
        function varargout = notify(O,varargin)
            varargout = notify@handle(O,varargin{:});
        end
        function varargout = listener(O,varargin)
            varargout = listener@handle(O,varargin{:});
        end
        function varargout = eq(O,varargin)
            varargout = eq@handle(O,varargin{:});
        end
        function varargout = ge(O,varargin)
            varargout = ge@handle(O,varargin{:});
        end
        function varargout = gt(O,varargin)
            varargout = gt@handle(O,varargin{:});
        end       
        function varargout = le(O,varargin)
            varargout = le@handle(O,varargin{:});
        end
        function varargout = lt(O,varargin)
            varargout = lt@handle(O,varargin{:});
        end
        function varargout = ne(O,varargin)
            varargout = ne@handle(O,varargin{:});
        end
    end
end
