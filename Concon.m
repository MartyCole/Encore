classdef Concon < matlab.mixin.Copyable
    properties (Access = private)        
        orig_st_points
        orig_en_points
        orig_st_idx
        orig_en_idx

        st_points
        en_points
        st_hemi
        en_hemi
        st_idx
        en_idx

        lh_P
        rh_P
        A  
        N

        lh_grid
        rh_grid
    end
    properties (Access = private, NonCopyable)
        lh_aabb
        rh_aabb
    end
    methods
        function obj = Concon(lh_grid, rh_grid, st_points, en_points, st_hemi, en_hemi)           
            % Grid theta, phi, and size P
            obj.lh_P = length(lh_grid.V);
            obj.rh_P = length(rh_grid.V);
            obj.A = [lh_grid.A; rh_grid.A] * [lh_grid.A; rh_grid.A].';

            obj.lh_grid = lh_grid;
            obj.rh_grid = rh_grid;
            obj.lh_aabb = AABBtree(lh_grid);
            obj.rh_aabb = AABBtree(rh_grid);

            obj.N = length(st_points);
            obj.st_hemi = st_hemi;
            obj.en_hemi = en_hemi;
            
            obj.st_points = zeros(obj.N,3);
            obj.en_points = zeros(obj.N,3);    
            obj.st_idx = zeros(obj.N,3);
            obj.en_idx = zeros(obj.N,3);   
            
            obj.get_coordinate_data(st_points, en_points)

            obj.orig_st_points = obj.st_points;
            obj.orig_en_points = obj.en_points;
            obj.orig_st_idx = obj.st_idx;
            obj.orig_en_idx = obj.en_idx;
        end

        function resample(obj, new_lh_grid, new_rh_grid)
            st_tmp_points = zeros(size(obj.st_points));

            tmp_idx = obj.st_idx(obj.st_hemi == 0,:);            
            st_tmp_points(obj.st_hemi == 0,:) = obj.st_points(obj.st_hemi == 0,1) .* obj.lh_grid.V(tmp_idx(:,1),:) + ...
                obj.st_points(obj.st_hemi == 0,2) .* obj.lh_grid.V(tmp_idx(:,2),:) + ...
                obj.st_points(obj.st_hemi == 0,3) .* obj.lh_grid.V(tmp_idx(:,3),:);

            tmp_idx = obj.st_idx(obj.st_hemi == 1,:) - obj.lh_P;
            st_tmp_points(obj.st_hemi == 1,:) = obj.st_points(obj.st_hemi == 1,1) .* obj.lh_grid.V(tmp_idx(:,1),:) + ...
                obj.st_points(obj.st_hemi == 1,2) .* obj.lh_grid.V(tmp_idx(:,2),:) + ...
                obj.st_points(obj.st_hemi == 1,3) .* obj.lh_grid.V(tmp_idx(:,3),:);

            en_tmp_points = zeros(size(obj.en_points));

            tmp_idx = obj.en_idx(obj.en_hemi == 0,:);            
            en_tmp_points(obj.en_hemi == 0,:) = obj.st_points(obj.en_hemi == 0,1) .* obj.rh_grid.V(tmp_idx(:,1),:) + ...
                obj.en_points(obj.en_hemi == 0,2) .* obj.rh_grid.V(tmp_idx(:,2),:) + ...
                obj.en_points(obj.en_hemi == 0,3) .* obj.rh_grid.V(tmp_idx(:,3),:);

            tmp_idx = obj.en_idx(obj.en_hemi == 1,:) - obj.lh_P;
            en_tmp_points(obj.en_hemi == 1,:) = obj.st_points(obj.en_hemi == 1,1) .* obj.rh_grid.V(tmp_idx(:,1),:) + ...
                obj.en_points(obj.en_hemi == 1,2) .* obj.rh_grid.V(tmp_idx(:,2),:) + ...
                obj.en_points(obj.en_hemi == 1,3) .* obj.rh_grid.V(tmp_idx(:,3),:);

            st_tmp_points = normr(st_tmp_points);
            en_tmp_points = normr(en_tmp_points);            

            obj.lh_P = length(new_lh_grid.V);
            obj.rh_P = length(new_rh_grid.V);
            obj.A = [new_lh_grid.A; new_rh_grid.A] * [new_lh_grid.A; new_rh_grid.A].';

            obj.lh_grid = new_lh_grid;
            obj.rh_grid = new_rh_grid;
            obj.lh_aabb = AABBtree(new_lh_grid);
            obj.rh_aabb = AABBtree(new_rh_grid);

            obj.get_coordinate_data(st_tmp_points, en_tmp_points);

            obj.orig_st_points = obj.st_points;
            obj.orig_en_points = obj.en_points;
            obj.orig_st_idx = obj.st_idx;
            obj.orig_en_idx = obj.en_idx;
        end

        function connectome = evaluate(obj, kernel, weighted)
            if weighted == true             
                [idx_a,idx_b] = ndgrid(1:3,1:3);

                % TODO: MAKE SURE THIS IS CORRECT FOR CREATING A WEIGHTED MATRIX                
                data = obj.st_points(:,idx_a(:)).*obj.en_points(:,idx_b(:));
                T_sp = obj.st_idx(:,idx_a(:));
                T_ep = obj.en_idx(:,idx_b(:));

                adjacency = accumarray([T_sp(:),T_ep(:)], data(:));
            else
                [~,idx_a] = min(obj.st_points,[],2,'linear');
                [~,idx_b] = min(obj.en_points,[],2,'linear');

                % create adjacency matrix from endpoints
                adjacency = sparse(obj.st_idx(idx_a), obj.en_idx(idx_b), 1, obj.lh_P + obj.rh_P, obj.lh_P + obj.rh_P);
            end

            adjacency = (adjacency + adjacency.') ./ (2*obj.N);

            % generate smooth connectome
            connectome = max(kernel * adjacency * kernel.',0);
        end

        function warp_connectome(obj,lh_warp,rh_warp)
            warp.V = [lh_warp.V; rh_warp.V];

            wx_sp = dot(obj.orig_st_points, reshape(warp.V(obj.orig_st_idx,1),obj.N,3),2);
            wy_sp = dot(obj.orig_st_points, reshape(warp.V(obj.orig_st_idx,2),obj.N,3),2);
            wz_sp = dot(obj.orig_st_points, reshape(warp.V(obj.orig_st_idx,3),obj.N,3),2);            
            
            warped_st_points = normr([wx_sp,wy_sp,wz_sp]);
            
            wx_ep = dot(obj.orig_en_points, reshape(warp.V(obj.orig_en_idx,1),obj.N,3),2);
            wy_ep = dot(obj.orig_en_points, reshape(warp.V(obj.orig_en_idx,2),obj.N,3),2);
            wz_ep = dot(obj.orig_en_points, reshape(warp.V(obj.orig_en_idx,3),obj.N,3),2);            
            
            warped_en_points = normr([wx_ep,wy_ep,wz_ep]);        

            obj.get_coordinate_data(warped_st_points, warped_en_points)
        end
    end

    methods(Access = protected)
        function cpObj = copyElement(obj)
            cpObj = copyElement@matlab.mixin.Copyable(obj);

            cpObj.lh_aabb = AABBtree(obj.lh_grid);
            cpObj.rh_aabb = AABBtree(obj.rh_grid);
        end
    end

    methods (Access = private)
        function get_coordinate_data(obj,st_points,en_points)                            
            [obj.st_points(obj.st_hemi == 0,:), ...
                obj.st_idx(obj.st_hemi == 0,:)] = obj.lh_aabb.get_barycentric_data(st_points(obj.st_hemi == 0,:),false);
            [obj.en_points(obj.en_hemi == 0,:), ...
                obj.en_idx(obj.en_hemi == 0,:)] = obj.lh_aabb.get_barycentric_data(en_points(obj.en_hemi == 0,:),false);
            [obj.st_points(obj.st_hemi == 1,:), ...
                obj.st_idx(obj.st_hemi == 1,:)] = obj.rh_aabb.get_barycentric_data(st_points(obj.st_hemi == 1,:),false);
            [obj.en_points(obj.en_hemi == 1,:), ...
                obj.en_idx(obj.en_hemi == 1,:)] = obj.rh_aabb.get_barycentric_data(en_points(obj.en_hemi == 1,:),false);
           
            obj.st_idx(obj.st_hemi == 1,:) = obj.st_idx(obj.st_hemi == 1,:) + obj.lh_P;
            obj.en_idx(obj.en_hemi == 1,:) = obj.en_idx(obj.en_hemi == 1,:) + obj.lh_P;
        end
    end
end

