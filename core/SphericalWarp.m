classdef SphericalWarp < matlab.mixin.Copyable
    % SphericalWarp
    %
    % This class represents a diffeomorphic warp (smooth, invertible mapping)
    % on the sphere. A warp is defined by a displacement field expressed in the
    % tangent plane at each vertex, and updated through:
    %
    %   - exponential map on the sphere
    %   - parallel transport of tangent vectors
    %   - barycentric interpolation using an AABB tree
    %
    % The class supports:
    %   - composing incremental warp updates
    %   - inverting a warp
    %   - applying rigid rotations
    %   - resampling a warp onto a new spherical grid
    %   - visualizing the warp field and its Jacobian determinant
    %
    % This class is used by Encore to apply and update diffeomorphic
    % transformations during spherical registration.
    %
    properties
        V % warped vertex positions
        T % triangulation
    end
    properties (Dependent)
        J % Jacobian determinant (area distortion)
    end
    properties (Access = private)        
        delta       
        P % number of vertices       
       
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
            %   - vertex positions V and triangulation T
            %   - original Voronoi areas for Jacobian computation
            %   - AABB tree for barycentric interpolation
            %
            % Input:
            %   grid : SphericalGrid object defining the base mesh

            obj.e1 = grid.e1;
            obj.e2 = grid.e2;
            obj.P = size(grid.V,1);

            % The warped vertices
            obj.V = grid.V;

            % The triangulation
            obj.T = grid.T;

            obj.grid_A = calc_voronoi_area(obj.V, obj.T);            
            obj.aabb = AABBtree(grid);
        end
        
        function obj = compose_warp(obj, displacement)
            % Compose an incremental displacement field with the current warp.
            %
            % Steps:
            %   - project displacement into the tangent plane (using e1, e2)
            %   - parallel transport displacement to the warped surface
            %   - interpolate displacement using barycentric coordinates
            %   - update vertex positions via the spherical exponential map
            %   - reject updates that would create folded triangles
            %
            % Input:
            %   displacement : Px2 tangent displacement field

            % Project displacement vector to the tangent plane in R3 then project to sphere 
            tangent_vec = (obj.e1 .* displacement(:,1) + obj.e2 .* displacement(:,2));            

            % Compose the new warp to the current overall warp
            [Vx,Tx] = obj.aabb.get_barycentric_data(obj.V, false);

            tangent_vec = obj.parallel_transport(tangent_vec(Tx,:)', obj.aabb.V(Tx,:)', repmat(obj.V,3,1)')';            
            tangent_vec = Vx(:) .* tangent_vec;
            tangent_vec = squeeze(sum(reshape(tangent_vec, obj.P, 3, 3), 2));
            
            % Project to sphere
            test_V = normr(sphere_exp_map(obj.V, tangent_vec));

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

	    function obj = invert_warp(obj)
            % Compute the inverse of the current warp.
            %
            % Uses:
            %   - logarithmic map to compute tangent vectors from warped to original grid
            %   - parallel transport to align tangent vectors
            %   - barycentric interpolation to reconstruct inverse displacement
            %   - exponential map to update vertex positions
            %
            % Produces a warp such that:  inverse_warp o warp ~ identity

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
            obj.V = normr(sphere_exp_map(obj.aabb.V, tangent_vec));
        end

        function obj = rotate(obj, rotation)
            % Apply a rigid rotation to the warp.
            %
            % Input:
            %   rotation : 3x3 rotation matrix (SO(3))
            %
            % Updates vertex positions as:
            %   V_new = V / rotation

            obj.V = obj.V / rotation;
        end

        function obj = resample(obj, new_grid)
            % Resample the warp onto a new spherical grid.
            %
            % Steps:
            %   - compute tangent displacement between original and warped vertices
            %   - interpolate displacement onto the new grid using barycentric lookup
            %   - parallel transport displacement to new grid vertices
            %   - update vertex positions via exponential map
            %   - update tangent frame, triangulation, and AABB tree
            %
            % Input:
            %   new_grid : SphericalGrid object defining the new mesh

            tangent_vec = sphere_log_map(obj.aabb.V,obj.V);

            % Compose the new warp to the current overall warp
            [Vx,Tx] = obj.aabb.get_barycentric_data(new_grid.V, false);

            tangent_vec = obj.parallel_transport(tangent_vec(Tx,:)', obj.aabb.V(Tx,:)', repmat(new_grid.V,3,1)')';            
            tangent_vec = Vx(:) .* tangent_vec;
            tangent_vec = squeeze(sum(reshape(tangent_vec, size(new_grid.V,1), 3, 3), 2));
            
            % Project to sphere
            obj.V = normr(sphere_exp_map(new_grid.V, tangent_vec));

            obj.e1 = new_grid.e1;
            obj.e2 = new_grid.e2;
            obj.P = size(new_grid.V,1);

            % The triangulation
            obj.T = new_grid.T;
            
            obj.aabb = AABBtree(new_grid);          
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

            tangent_vec = sphere_log_map(obj.V,obj.aabb.V);                      

            quiver3(obj.aabb.V(:,1),obj.aabb.V(:,2),obj.aabb.V(:,3), ...
                tangent_vec(:,1), tangent_vec(:,2), tangent_vec(:,3),2, ...
                'LineWidth',2);
            hold on
            trisurf(obj.T, obj.V(:,1), obj.V(:,2), obj.V(:,3), 'FaceColor',[0.8,0.8,0.8], ...
                'EdgeAlpha', 0)
            hold off
            title(fig_title)
            axis off
            axis equal
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
            pq = dot(orig_pt, new_pt, 1);
            vq = dot(tan_vec, new_pt, 1);

            denom = 1 + pq;

            % compute transported vector
            correction = (vq ./ denom) .* (orig_pt + new_pt);
            new_vec = tan_vec - correction;

            % handle antipodal or near-antipodal cases
            idx_bad = denom < 1e-8;
            new_vec(:, idx_bad) = tan_vec(:, idx_bad);
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
