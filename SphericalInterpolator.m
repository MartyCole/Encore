classdef SphericalInterpolator < handle
    properties (SetAccess = immutable, GetAccess = private)
        tree_ptr uint64;
        V
        T
    end
    properties (Access = private)          
        Vq 
        Tq    
        initialised
    end
    methods
        function obj = SphericalInterpolator(grid, varargin)            
            obj.V = grid.V;
            obj.T = grid.T - 1;

            % initialise the AABB tree
            obj.tree_ptr = init_AABB_tree_mex(obj.V, obj.T);

            if nargin == 2                
                obj.update_query_points(varargin{1});   
                obj.initialised = true;
            else
                obj.initialised = false;
            end            
        end

        function [new_Q] = evaluate_2d(obj, Q)
            if obj.initialised == false
                error('Interpolation object not yet initialised')
            end

            if isreal(Q)
                new_Q = bary_interp_2D_mex(obj.Vq,obj.Vq,obj.Tq,obj.Tq,Q);            
            else
                error('Cannot interpolate complex values!')
            end
        end   

        function update_query_points(obj, points)
            % get barycentric coordinates 
            [obj.Vq,obj.Tq] = query_AABB_tree_mex(obj.tree_ptr, obj.V, obj.T, points);
           
            % save coordinates for interpolating data
            obj.Vq = obj.Vq.'; 
            obj.Tq = int64(obj.Tq).';            

            obj.initialised = true;
        end

        function update_query_points_inv(obj, points)
            % get barycentric coordinates 
            [obj.Vq,obj.Tq] = query_AABB_tree_mex(obj.tree_ptr, points, obj.T, obj.V);
           
            % save coordinates for interpolating data
            obj.Vq = obj.Vq.'; 
            obj.Tq = int64(obj.Tq).';            

            obj.initialised = true;
        end

        function [Vq,Tq] = get_barycentric_data(obj, points)
            % get barycentric coordinates 
            [Vq,Tq] = query_AABB_tree_mex(obj.tree_ptr, obj.V, obj.T, points); 
            Tq = int64(Tq) + 1;
        end

        function delete(obj)
            destroy_AABB_tree_mex(obj.tree_ptr);
        end
    end    
end

