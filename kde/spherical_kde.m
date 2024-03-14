function connectome = spherical_kde(grid,kernel,pts_in,pts_out,hemi_in,hemi_out,weighted)

% find closet vertex for each point
tree = AABBtree(grid);
[D_sp,T_sp] = tree.get_barycentric_data(pts_in,false);
[D_ep,T_ep] = tree.get_barycentric_data(pts_out,false);

T_sp = T_sp + int64(hemi_in.' * size(grid.V, 1));
T_ep = T_ep + int64(hemi_out.' * size(grid.V, 1));

if weighted == true             
    [idx_a,idx_b] = ndgrid(1:3,1:3);

    data = D_sp(:,idx_a(:)).*D_ep(:,idx_b(:));
    T_sp = T_sp(:,idx_a(:));
    T_ep = T_ep(:,idx_b(:));

    A = accumarray([T_sp(:),T_ep(:)], data(:));
else
    [~,idx_a] = min(D_sp,[],2,'linear');
    [~,idx_b] = min(D_ep,[],2,'linear');
    
    % create adjacency matrix from endpoints
    A = full(sparse(T_sp(idx_a), T_ep(idx_b), 1, 2*size(grid.V,1), 2*size(grid.V,1)));
end

A = A + A.';

% generate smooth connectome
connectome = max(kernel * A * kernel.',0);

end

