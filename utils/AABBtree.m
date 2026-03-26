classdef AABBtree < handle
    properties (SetAccess = private, GetAccess = private)
        tree_ptr
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
            obj.tree_ptr = aabb_mex('init', obj.V, obj.T);
        end     

        function [Vq,Tq] = get_barycentric_data(obj, points, CPPformat)
            if isempty(obj.tree_ptr) 
                error('AABBTree: tree has been destroyed'); 
            end

            % get barycentric coordinates 
            [Vq,Tq] = aabb_mex('query', obj.tree_ptr, points); 

            if CPPformat
                Vq = Vq.'; 
                Tq = int64(Tq).'; 
            else
                Tq = int64(Tq) + 1;
            end
        end

        function delete(obj)
            % Destructor: free the tree 
            if ~isempty(obj.tree_ptr) 
                try 
                    aabb_mex('destroy', obj.tree_ptr); 
                catch ME 
                    warning(ME.identifier, 'AABBTree: destroy failed: %s', ME.message); 
                end 
                
                obj.tree_ptr = []; 
            end
        end
    end
end

