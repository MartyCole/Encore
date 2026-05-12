classdef SConcon < Concon
    % SConcon    Continuous structural connectivity representation
    %
    % This class represents a structural connectome as pairs of endpoints
    % lying on two unit icosahedral sphere meshes (left and right hemispheres).
    % Each fiber endpoint is expressed in barycentric coordinates relative to
    % a triangle on the mesh. The class supports:
    %   - computing adjacency matrices
    %   - evaluating smooth connectomes using kernels
    %   - warping endpoints to new spherical meshes
    %   - rotating matrices using barycentric interpolation
    %
    % The class stores both original and warped endpoint coordinates and
    % triangle indices, enabling repeated warping without loss of precision.  
    %
    properties (Access = private)        
        orig_st_points % Original (unwarped) barycentric data
        orig_en_points % Original (unwarped) barycentric data
        orig_st_idx    % Original (unwarped) barycentric data
        orig_en_idx    % Original (unwarped) barycentric data   

        N
    end

    properties (SetAccess = private)        
        st_points      % Nx3 barycentric coordinates of fiber start points
        en_points      % Nx3 barycentric coordinates of fiber end points
        st_hemi        % Nx1 hemisphere indicator for start points (0=LH, 1=RH)
        en_hemi        % Nx1 hemisphere indicator for end points (0=LH, 1=RH)
        st_idx         % Nx3 triangle vertex indices for start points
        en_idx         % Nx3 triangle vertex indices for end points
    end
 
    methods
        function obj = SConcon(lh_grid, rh_grid, st_points, en_points, st_hemi, en_hemi, gpu)     
            % Constructor for SConcon
            %
            % Inputs:
            %   lh_grid, rh_grid : spherical meshes for LH and RH
            %   st_points, en_points : Nx3 Cartesian coordinates of fiber endpoints
            %   st_hemi, en_hemi : hemisphere flags for each endpoint
            %   gpu: (true) enable gpu acceleration

            if nargin < 7
                gpu = true;
            end
            
            % initialise this Concon object using the given grid and endpoints
            obj.initialize(lh_grid, rh_grid, st_points, en_points, st_hemi, en_hemi, gpu);
        end

        function obj = reset(obj)
            % Reset endpoints to last applied warp

            obj.st_points = obj.orig_st_points;
            obj.en_points = obj.orig_en_points;
            obj.st_idx = obj.orig_st_idx;
            obj.en_idx = obj.orig_en_idx;
        end     

        function adjacency = get_adjacency(obj)
            % Compute the fiber adjacency matrix.
            %
            % Uses a compiled MEX function for speed. The adjacency matrix counts
            % contributions of fibers between mesh triangles, normalized by N.
            %
            % Output is returned as a symmetric gpuArray.

            P = obj.lh_P + obj.rh_P;

            % calculate adjacency matrix (must be done on CPU) using mex for speed
            adjacency = build_adjacency(int32(obj.st_idx), int32(obj.en_idx), ...
                                            obj.st_points, obj.en_points, P, P);    

            % force symmetry and normalise (can move back to the GPU now) 
            adjacency = obj.maybe_GPU((adjacency + adjacency.') ./ (2*obj.N));
        end

        function [connectome, De1, De2] = evaluate(obj, kernel, d_kernel)
            % Evaluate the smooth connectome using a kernel defined on the sphere.
            %
            % connectome = K' * A * K   (clamped to nonnegative)
            %
            % If derivative kernels are provided, compute directional derivatives
            % projected onto the tangent basis (e1, e2) of each vertex. Forced
            % to always non-negative, ensure Kernel supports this.

            if nargin < 3
                [connectome, De1, De2] = evaluate@Concon(obj, kernel);
            else
                [connectome, De1, De2] = evaluate@Concon(obj, kernel, d_kernel);           
            end         

            % force the connectome to be non-negative
            connectome = max(connectome,0);
        end

        function [Q, Q_e1, Q_e2] = Q_transform(obj, kernel, d_kernel)
            % Convert F to Q = sqrt(F) and propagate derivatives via chain rule.

            if nargin > 2
                % evaluate F and its derivatives
                [F, F_e1, F_e2] = obj.evaluate(kernel, d_kernel);  

                % perform the transformation
                Q = sqrt(F);
                
                % propegate the square root to the derivative by chain rule
                Q_e1 = ((1./(2*(max(Q,1e-15)))) .* F_e1);
                Q_e2 = ((1./(2*(max(Q,1e-15)))) .* F_e2);
            else
                % evaluate only F
                Q = sqrt(obj.evaluate(kernel)); 
            end            
        end
        
        function [warped_st_points,warped_en_points] = warp_connectome(obj, lh_warp, rh_warp)
            % Warp fiber endpoints to a new spherical mesh.
            %
            % Uses barycentric coordinates to interpolate new positions, then
            % reprojects onto the unit sphere and recomputes barycentric data.

            warp.V = [lh_warp.V; rh_warp.V];

            % warp the start points by using the barycentric coordinates on the new grid
            wx_sp = dot(obj.orig_st_points, reshape(warp.V(obj.orig_st_idx,1),obj.N,3),2);
            wy_sp = dot(obj.orig_st_points, reshape(warp.V(obj.orig_st_idx,2),obj.N,3),2);
            wz_sp = dot(obj.orig_st_points, reshape(warp.V(obj.orig_st_idx,3),obj.N,3),2);            
            
            % warp the end points by using the barycentric coordinates on the new grid
            wx_ep = dot(obj.orig_en_points, reshape(warp.V(obj.orig_en_idx,1),obj.N,3),2);
            wy_ep = dot(obj.orig_en_points, reshape(warp.V(obj.orig_en_idx,2),obj.N,3),2);
            wz_ep = dot(obj.orig_en_points, reshape(warp.V(obj.orig_en_idx,3),obj.N,3),2);            
            
            % make sure the points are on the unit sphere
            warped_st_points = normr([wx_sp,wy_sp,wz_sp]);
            warped_en_points = normr([wx_ep,wy_ep,wz_ep]);        

            % calculate the new barycentric coordinates in relation to the original grid
            obj.get_coordinate_data(warped_st_points, warped_en_points)
        end        

        function obj = apply_warp(obj, new_lh_grid, new_rh_grid)
            % Apply a warp to the current mesh or reinitialize on a new mesh.
            %
            % If no new grids are provided, the warp is applied in-place.
            % Otherwise, the Concon is fully reinitialized using the warped points.

            warp.V = [obj.lh_grid.V; obj.rh_grid.V];

            % warp the start points by using the barycentric coordinates on the new grid
            wx_sp = dot(obj.st_points, reshape(warp.V(obj.st_idx,1),obj.N,3),2);
            wy_sp = dot(obj.st_points, reshape(warp.V(obj.st_idx,2),obj.N,3),2);
            wz_sp = dot(obj.st_points, reshape(warp.V(obj.st_idx,3),obj.N,3),2);            
            
            % warp the end points by using the barycentric coordinates on the new grid
            wx_ep = dot(obj.en_points, reshape(warp.V(obj.en_idx,1),obj.N,3),2);
            wy_ep = dot(obj.en_points, reshape(warp.V(obj.en_idx,2),obj.N,3),2);
            wz_ep = dot(obj.en_points, reshape(warp.V(obj.en_idx,3),obj.N,3),2);            
            
            % make sure the points are on the unit sphere
            warped_st_points = normr([wx_sp,wy_sp,wz_sp]);
            warped_en_points = normr([wx_ep,wy_ep,wz_ep]);        

            if nargin == 1
                % apply warp to original grid if no new ones given
                obj.get_coordinate_data(warped_st_points, warped_en_points)

                obj.orig_st_points = obj.st_points;
                obj.orig_en_points = obj.en_points;
                obj.orig_st_idx = obj.st_idx;
                obj.orig_en_idx = obj.en_idx;
            else
                % re-initialise the current Concon object using the new grid
                obj.initialize(new_lh_grid, new_rh_grid, warped_st_points, warped_en_points, obj.st_hemi, obj.en_hemi)
            end
        end

    end

    methods (Access = protected)        
        function initialize(obj, lh_grid, rh_grid, st_points, en_points, st_hemi, en_hemi, gpu)   
            % Core initialization routine.

            initialize@Concon(obj, lh_grid, rh_grid, gpu);

            % endpoint information
            obj.N = length(st_points);
            obj.st_hemi = int8(st_hemi);
            obj.en_hemi = int8(en_hemi);
            
            % empty matrices to fill with barycentric information (get_coordinate_data)
            obj.st_points = zeros(obj.N,3);
            obj.en_points = zeros(obj.N,3);    
            obj.st_idx = zeros(obj.N,3,'int16');
            obj.en_idx = zeros(obj.N,3,'int16');   
            
            % populate the barycentric matrices
            obj.get_coordinate_data(st_points, en_points)

            % save endpoint data before warping
            obj.orig_st_points = obj.st_points;
            obj.orig_en_points = obj.en_points;
            obj.orig_st_idx = obj.st_idx;
            obj.orig_en_idx = obj.en_idx;
        end

        function get_coordinate_data(obj, st_points, en_points)  
            % Compute barycentric coordinates and triangle indices for endpoints.

            % create masks for each hemisphere and start/end point set
            lh_mask_st = (obj.st_hemi == 0); rh_mask_st = ~lh_mask_st;
            lh_mask_en = (obj.en_hemi == 0); rh_mask_en = ~lh_mask_en;
            
            % extract the subsets of points to query on each hemisphere grid
            st_lh = st_points(lh_mask_st, :); st_rh = st_points(rh_mask_st, :);
            en_lh = en_points(lh_mask_en, :); en_rh = en_points(rh_mask_en, :);
            
            % calculate the barycentric coordinates of the given points on the original grid
            [st_lh_out, st_lh_idx] = obj.lh_aabb.get_barycentric_data(st_lh, false);
            [en_lh_out, en_lh_idx] = obj.lh_aabb.get_barycentric_data(en_lh, false);
            
            % calculate the barycentric coordinates of the given points on the original grid
            [st_rh_out, st_rh_idx] = obj.rh_aabb.get_barycentric_data(st_rh, false);
            [en_rh_out, en_rh_idx] = obj.rh_aabb.get_barycentric_data(en_rh, false);
            
            % write the subsets barycentric coordinates back to the original data objects
            obj.st_points(lh_mask_st, :) = st_lh_out; obj.st_idx(lh_mask_st, :) = st_lh_idx;
            obj.en_points(lh_mask_en, :) = en_lh_out; obj.en_idx(lh_mask_en, :) = en_lh_idx;
            
            % the right hemisphere starts at a higher index, so we add that here
            obj.st_points(rh_mask_st, :) = st_rh_out; obj.st_idx(rh_mask_st, :) = st_rh_idx + obj.lh_P;
            obj.en_points(rh_mask_en, :) = en_rh_out; obj.en_idx(rh_mask_en, :) = en_rh_idx + obj.lh_P;            
        end

    end

end

