function [warped_sp,warped_ep] = warp_connectome(grid,warp,sp,ep)

tree = AABBtree(grid);
[D_sp,T_sp] = tree.get_barycentric_data(sp',false);
[D_ep,T_ep] = tree.get_barycentric_data(ep',false);

wx_sp = dot(D_sp, reshape(warp.V(T_sp,1),length(sp),3),2);
wy_sp = dot(D_sp, reshape(warp.V(T_sp,2),length(sp),3),2);
wz_sp = dot(D_sp, reshape(warp.V(T_sp,3),length(sp),3),2);            

warped_sp = normr([wx_sp,wy_sp,wz_sp]).';

wx_ep = dot(D_ep, reshape(warp.V(T_ep,1),length(ep),3),2);
wy_ep = dot(D_ep, reshape(warp.V(T_ep,2),length(ep),3),2);
wz_ep = dot(D_ep, reshape(warp.V(T_ep,3),length(ep),3),2);            

warped_ep = normr([wx_ep,wy_ep,wz_ep]).';

end