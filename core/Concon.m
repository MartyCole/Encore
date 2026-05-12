classdef Concon < matlab.mixin.Copyable
    % Concon    Continuous connectivity representation
    %
    % This class represents the Abstract class to support:
    %   - computing adjacency matrices
    %   - evaluating smooth connectomes using kernels
    %   - warping endpoints to new spherical meshes
    %   - rotating matrices using barycentric interpolation
    %
    % The class stores both original and warped endpoint coordinates and
    % triangle indices, enabling repeated warping without loss of precision.  
    %
    properties
        use_GPU logical = true % convert matrices to GPU when enabled
    end
    properties (Access = protected)      
        lh_P           % Number of vertices in LH/RH meshes (lh_P + rh_P = P)
        rh_P           % Number of vertices in LH/RH meshes (lh_P + rh_P = P)        
        
        lh_grid        % Left hemisphere spherical mesh
        rh_grid        % Right hemisphere spherical mesh
    end
    properties (SetAccess = protected)        
        A              % PxP matrix of triangle areas used for integration
    end
    properties (Access = protected, NonCopyable)
        lh_aabb        % AABB tree for fast LH point‑location queries
        rh_aabb        % AABB tree for fast RH point‑location queries
    end

    methods (Abstract)    
        % Reset endpoints to last applied warp
        obj = reset(obj)

        % Compute the ConCon adjacency matrix.
        adjacency = get_adjacency(obj)

        % Apply a warp to the current mesh or reinitialize on a new mesh.
        obj = apply_warp(obj, new_lh_grid, new_rh_grid)
        
        % Warp ConCon to a new spherical mesh.
        [warped_st_points,warped_en_points] = warp_connectome(obj, lh_warp, rh_warp)

        % Transformation of the Concon to Hilbert Sphere S_inf
        [Q, Q_e1, Q_e2] = Q_transform(obj, kernel, kernel_diff)
    end

    methods
        function [connectome, De1, De2] = evaluate(obj, kernel, d_kernel)  
            % Evaluate the smooth connectome using a kernel defined on the sphere.
            %
            % connectome = K' * A * K   (clamped to nonnegative)
            %
            % If derivative kernels are provided, compute directional derivatives
            % projected onto the tangent basis (e1, e2) of each vertex.

            % get the adjacency matrix
            adjacency = obj.get_adjacency();

            % generate the smooth connectome using the given kernel
            connectome = kernel.' * adjacency * kernel;

            De1 = [];
            De2 = [];

            % generate the derivative (Dk * A * K' + K * A * Dk') if needed
            if nargin == 3                        
                AKt = adjacency * kernel.';
    
                Dx = (d_kernel.x * AKt); Dx = Dx + Dx.';
                Dy = (d_kernel.y * AKt); Dy = Dy + Dy.';
                Dz = (d_kernel.z * AKt); Dz = Dz + Dz.';
    
                % project to the tangent space of theta and phi
                e1 = single([obj.lh_grid.e1; obj.rh_grid.e1]);
                e2 = single([obj.lh_grid.e2; obj.rh_grid.e2]);
                De1 = Dx.*e1(:,1) + Dy.*e1(:,2) + Dz.*e1(:,3);
                De2 = Dx.*e2(:,1) + Dy.*e2(:,2) + Dz.*e2(:,3);
            end
        end

        function F_rot = rotate_matrix_barycentric(obj, F, R_lh, R_rh)
            % Rotate a matrix defined on the spherical mesh using barycentric interpolation.

            % rotate vertices
            V_lh_rot = (R_lh * obj.lh_grid.V.').';
            V_rh_rot = (R_rh * obj.rh_grid.V.').';
        
            % barycentric data on original grids
            [bc_lh, tri_lh] = obj.lh_aabb.get_barycentric_data(V_lh_rot, false);
            [bc_rh, tri_rh] = obj.rh_aabb.get_barycentric_data(V_rh_rot, false);
        
            % left hemisphere
            rows_lh = repmat((1:obj.lh_P).', 3,1);      
            cols_lh = tri_lh(:);                  
            vals_lh = bc_lh(:);                    
            W_lh = sparse(rows_lh, cols_lh, vals_lh, obj.lh_P, obj.lh_P);
        
            % right hemisphere
            rows_rh = repmat((1:obj.rh_P).', 3,1);
            cols_rh = tri_rh(:);
            vals_rh = bc_rh(:);
            W_rh = sparse(rows_rh, cols_rh, vals_rh, obj.rh_P, obj.rh_P);        
         
            % block-diagonal interpolation
            W = obj.maybe_GPU((blkdiag(W_lh, W_rh)));
        
            % rotate full matrix
            F_rot = single(W * double(F) * W.');
        end
    end

    methods(Access = protected)
        function initialize(obj, lh_grid, rh_grid, gpu)   
            % Core initialization routine.

            obj.use_GPU = gpu;

            % the grids to evaluate from
            obj.lh_grid = lh_grid;
            obj.rh_grid = rh_grid;

            % the size of the grids
            obj.lh_P = int64(length(lh_grid.V)); 
            obj.rh_P = int64(length(rh_grid.V));

            % area of triangles for use in integration
            obj.A = [lh_grid.A; rh_grid.A] * [lh_grid.A; rh_grid.A].';

            % AABB trees for each grid
            obj.lh_aabb = AABBtree(lh_grid);
            obj.rh_aabb = AABBtree(rh_grid);
        end

        function cpObj = copyElement(obj)
            % Deep copy method for Copyable mixin.
            %
            % Ensures AABB trees are recreated rather than shallow‑copied.

            cpObj = copyElement@matlab.mixin.Copyable(obj);

            cpObj.lh_aabb = AABBtree(obj.lh_grid);
            cpObj.rh_aabb = AABBtree(obj.rh_grid);
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

