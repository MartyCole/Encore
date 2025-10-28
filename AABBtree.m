classdef AABBtree < handle
    properties (SetAccess = private, GetAccess = private)
        tree_ptr uint64;
    end

    properties (SetAccess = private, GetAccess = public)
        V
        T
    end

    methods
        function obj = AABBtree(grid)            
            obj.V = grid.V;
            obj.T = grid.T - 1;

            % initialise the AABB tree
            obj.tree_ptr = init_AABB_tree_mex(obj.V, obj.T);
        end     

        function [Vq,Tq] = get_barycentric_data(obj, points, CPPformat)
            % get barycentric coordinates 
            [Vq,Tq] = query_AABB_tree_mex(obj.tree_ptr, obj.V, obj.T, points); 

            if CPPformat
                Vq = Vq.'; 
                Tq = int64(Tq).'; 
            else
                Tq = int64(Tq) + 1;
            end
        end

        function delete(obj)
            destroy_AABB_tree_mex(obj.tree_ptr);
        end
    end

    methods (Static = true)
        function obj = loadobj(obj)
            % initialise the AABB tree
            obj.tree_ptr = init_AABB_tree_mex(obj.V, obj.T);
        end
    end    
end

