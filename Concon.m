classdef Concon
    properties (Access = private)        
        delta
        
        P  
        P2
        A        

        lh_aabb
        rh_aabb

        D_sp
        D_ep
        T_sp
        T_ep
        in_idx
        out_idx

        lh_dst
        rh_dst
        kernel
    end
    methods
        function obj = Concon(lh_grid, rh_grid, pts_in, pts_out, FWHM, epsilon, delta)          
            % Small value for derivatives
            obj.delta = delta*2;   
            
            % Grid theta, phi, and size P
            obj.P = length(lh_grid.V);
            obj.P2 = length(lh_grid.V) + length(rh_grid.V);
            obj.A = [lh_grid.A; rh_grid.A] * [lh_grid.A; rh_grid.A].';

            obj.lh_aabb = AABBtree(lh_grid);
            obj.rh_aabb = AABBtree(rh_grid);               

            % Barycentric coords for endpoints
            obj.in_idx = pts_in(:,4) == 0;
            obj.out_idx = pts_out(:,4) == 0;

            [obj.D_sp(obj.in_idx,:),obj.T_sp(obj.in_idx,:)] = obj.lh_aabb.get_barycentric_data(pts_in(obj.in_idx,1:3),false);
            [obj.D_ep(obj.out_idx,:),obj.T_ep(obj.out_idx,:)] = obj.lh_aabb.get_barycentric_data(pts_out(obj.out_idx,1:3),false);
            [obj.D_sp(~obj.in_idx,:),obj.T_sp(~obj.in_idx,:)] = obj.rh_aabb.get_barycentric_data(pts_in(~obj.in_idx,1:3),false);
            [obj.D_ep(~obj.out_idx,:),obj.T_ep(~obj.out_idx,:)] = obj.rh_aabb.get_barycentric_data(pts_out(~obj.out_idx,1:3),false); 

            % radius of an idealised spherical brain
            R = 180/pi;    
            
            % get geodesic distance between each pair of vertices on the grid
            obj.lh_dst = acos(max(min(lh_grid.V * lh_grid.V.',1),-1))*R;
            obj.rh_dst = acos(max(min(rh_grid.V * rh_grid.V.',1),-1))*R;
            
            obj.kernel = obj.spherical_kernel(FWHM,epsilon);
        end    

        function F = evaluate(obj, varargin)
            % find closet vertex for each point  
            new_D_sp = obj.D_sp;
            new_D_ep = obj.D_ep;
            new_T_sp = obj.T_sp;
            new_T_ep = obj.T_ep;
            
            if nargin == 3
                [pts_in, pts_out] = obj.warp(varargin{1},varargin{2});

                [new_D_sp(obj.in_idx,:),new_T_sp(obj.in_idx,:)] = obj.lh_aabb.get_barycentric_data(pts_in(obj.in_idx,:),false);
                [new_D_ep(obj.out_idx,:),new_T_ep(obj.out_idx,:)] = obj.lh_aabb.get_barycentric_data(pts_out(obj.out_idx,:),false);
            
                [new_D_sp(~obj.in_idx,:),new_T_sp(~obj.in_idx,:)] = obj.rh_aabb.get_barycentric_data(pts_in(~obj.in_idx,:),false);
                [new_D_ep(~obj.out_idx,:),new_T_ep(~obj.out_idx,:)] = obj.rh_aabb.get_barycentric_data(pts_out(~obj.out_idx,:),false);
            end
            
            new_T_sp(~obj.in_idx,:) = new_T_sp(~obj.in_idx,:) + obj.P;
            new_T_ep(~obj.out_idx,:) = new_T_ep(~obj.out_idx,:) + obj.P;
            
            % create the weighted adjacency matrix
            adj = accumarray([new_T_sp(:),new_T_ep(:)],new_D_sp(:),[obj.P2,obj.P2]) + ...
                    accumarray([new_T_ep(:),new_T_sp(:)],new_D_ep(:),[obj.P2,obj.P2]);
            
            % generate smooth connectome
            F = obj.kernel * adj * obj.kernel';
            F(abs(F) < 1e-16) = 0;
            F(F < 0) = 0;
            F = (F ./ sum(F(:) .* obj.A(:)));
        end
    end

    methods (Access = private)
        function [pts_in, pts_out] = warp(obj, lh_warp, rh_warp)
            wx_sp = zeros(length(obj.in_idx),1);
            wy_sp = zeros(length(obj.in_idx),1);
            wz_sp = zeros(length(obj.in_idx),1);
            wx_ep = zeros(length(obj.out_idx),1);
            wy_ep = zeros(length(obj.out_idx),1);
            wz_ep = zeros(length(obj.out_idx),1);
                                   
            wx_sp(obj.in_idx) = dot(obj.D_sp(obj.in_idx,:), reshape(lh_warp.V(obj.T_sp(obj.in_idx,:),1),[],3),2);
            wy_sp(obj.in_idx) = dot(obj.D_sp(obj.in_idx,:), reshape(lh_warp.V(obj.T_sp(obj.in_idx,:),2),[],3),2);
            wz_sp(obj.in_idx) = dot(obj.D_sp(obj.in_idx,:), reshape(lh_warp.V(obj.T_sp(obj.in_idx,:),3),[],3),2);     
            
            wx_ep(obj.out_idx) = dot(obj.D_ep(obj.out_idx,:), reshape(lh_warp.V(obj.T_ep(obj.out_idx,:),1),[],3),2);
            wy_ep(obj.out_idx) = dot(obj.D_ep(obj.out_idx,:), reshape(lh_warp.V(obj.T_ep(obj.out_idx,:),2),[],3),2);
            wz_ep(obj.out_idx) = dot(obj.D_ep(obj.out_idx,:), reshape(lh_warp.V(obj.T_ep(obj.out_idx,:),3),[],3),2);            

            wx_sp(~obj.in_idx) = dot(obj.D_sp(~obj.in_idx,:), reshape(rh_warp.V(obj.T_sp(~obj.in_idx,:),1),[],3),2);
            wy_sp(~obj.in_idx) = dot(obj.D_sp(~obj.in_idx,:), reshape(rh_warp.V(obj.T_sp(~obj.in_idx,:),2),[],3),2);
            wz_sp(~obj.in_idx) = dot(obj.D_sp(~obj.in_idx,:), reshape(rh_warp.V(obj.T_sp(~obj.in_idx,:),3),[],3),2);     
            
            wx_ep(~obj.out_idx) = dot(obj.D_ep(~obj.out_idx,:), reshape(rh_warp.V(obj.T_ep(~obj.out_idx,:),1),[],3),2);
            wy_ep(~obj.out_idx) = dot(obj.D_ep(~obj.out_idx,:), reshape(rh_warp.V(obj.T_ep(~obj.out_idx,:),2),[],3),2);
            wz_ep(~obj.out_idx) = dot(obj.D_ep(~obj.out_idx,:), reshape(rh_warp.V(obj.T_ep(~obj.out_idx,:),3),[],3),2);

            pts_in = normr([wx_sp,wy_sp,wz_sp]);    
            pts_out = normr([wx_ep,wy_ep,wz_ep]);
        end

        function kernel = spherical_kernel(obj,FWHM,epsilon)   
            kernel = zeros(obj.P2,obj.P2);

            % calculate the bandwidth and truncation distance
            sigma = (FWHM/2) / sqrt(8*log(2));
            max_dist = (FWHM/2) * sqrt((-log2(epsilon)));
                        
            % truncate anything greater than max_dist
            obj.lh_dst(obj.lh_dst > max_dist) = Inf;
            obj.rh_dst(obj.rh_dst > max_dist) = Inf;

            % calculate the kernel for the grid
            kernel(1:obj.P,1:obj.P) = exp(-(obj.lh_dst.^2 / (2 * (sigma^2))));
            kernel((obj.P+1):end,(obj.P+1):end) = exp(-(obj.rh_dst.^2 / (2 * (sigma^2))));

            % normalise the kernel so that counts stay the same
            kernel = sparse(kernel ./ sum(kernel,2));
        end
    end
end

