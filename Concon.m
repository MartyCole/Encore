classdef Concon
    properties (Access = private)        
        delta
        
        P    
        A
        
        p_Vx
        p_Tx
        p_Vy
        p_Ty
        m_Vx
        m_Tx
        m_Vy
        m_Ty
        
        Vg
        Tg

        lh_aabb
        rh_aabb
    end
    methods
        function obj = Concon(lh_grid, rh_grid, delta)          
            % Small value for derivatives
            obj.delta = delta*2;   
            
            % Grid theta, phi, and size P
            obj.P = length(lh_grid.V);
            obj.A = [lh_grid.A; rh_grid.A] * [lh_grid.A; rh_grid.A].';

            obj.lh_aabb = AABBtree(lh_grid);
            obj.rh_aabb = AABBtree(rh_grid);
            
            % Interpolated values at delta coordinates
            [obj.p_Vx,obj.p_Tx] = obj.get_coordinate_data(sphere_exp_map(lh_grid.V, lh_grid.e1 * delta), ...
                                                          sphere_exp_map(rh_grid.V, rh_grid.e1 * delta));
            [obj.m_Vx,obj.m_Tx] = obj.get_coordinate_data(sphere_exp_map(lh_grid.V, lh_grid.e1 * -delta), ...
                                                          sphere_exp_map(rh_grid.V, rh_grid.e1 * -delta));
            [obj.p_Vy,obj.p_Ty] = obj.get_coordinate_data(sphere_exp_map(lh_grid.V, lh_grid.e2 * delta), ...
                                                          sphere_exp_map(rh_grid.V, rh_grid.e2 * delta));
            [obj.m_Vy,obj.m_Ty] = obj.get_coordinate_data(sphere_exp_map(lh_grid.V, lh_grid.e2 * -delta), ...
                                                          sphere_exp_map(rh_grid.V, rh_grid.e2 * -delta));  
            
            [obj.Vg,obj.Tg] = obj.get_coordinate_data(lh_grid.V,rh_grid.V);            
        end

        function [dQe1,dQe2] = get_derivative(obj, Q)          
            dQe1 = (bary_interp_2D_mex(obj.Vg,obj.p_Vx,obj.Tg,obj.p_Tx,Q) - ...
                    bary_interp_2D_mex(obj.Vg,obj.m_Vx,obj.Tg,obj.m_Tx,Q)) / obj.delta;
            dQe2 = (bary_interp_2D_mex(obj.Vg,obj.p_Vy,obj.Tg,obj.p_Ty,Q) - ...
                    bary_interp_2D_mex(obj.Vg,obj.m_Vy,obj.Tg,obj.m_Ty,Q)) / obj.delta;
        end     

        function new_Q = evaluate_Q(obj, Q, lh_warp, rh_warp)
            dJ = [lh_warp.J;rh_warp.J];
            dJ = sqrt(dJ) * sqrt(dJ).';     

            [Vq,Tq] = obj.get_coordinate_data(lh_warp.V, rh_warp.V);  
            new_Q = bary_interp_2D_mex(Vq,Vq,Tq,Tq,Q) .* dJ;     
            new_Q = (new_Q + new_Q.') / 2;                
            new_Q = new_Q / sqrt(sum(new_Q(:).^2 .* obj.A(:)));
        end      

        function new_F = evaluate(obj, F, lh_warp, rh_warp)
            dJ = [lh_warp.J;rh_warp.J];
            dJ = dJ * dJ.';     

            [Vq,Tq] = obj.get_coordinate_data(lh_warp.V, rh_warp.V);  
            new_F = bary_interp_2D_mex(Vq,Vq,Tq,Tq,F) .* dJ;     
            new_F = (new_F + new_F.') / 2;                
        end      
    end

    methods (Access = private)
        function [dV,dT] = get_coordinate_data(obj,lh_pts,rh_pts)                            
            % Interpolated values at delta coordinates
            [lh_V,lh_T] = obj.lh_aabb.get_barycentric_data(lh_pts, true);
            [rh_V,rh_T] = obj.rh_aabb.get_barycentric_data(rh_pts, true);
            
            dV = [lh_V,rh_V];
            dT = [lh_T,rh_T + obj.P];
        end
    end
end

