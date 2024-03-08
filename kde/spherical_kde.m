function connectome = spherical_kde(grid,kernel,pts_in,pts_out,hemi_in,hemi_out)

% find closet vertex for each point
tree = AABBtree(grid);
[D_sp,T_sp] = tree.get_barycentric_data(pts_in,false);
[D_ep,T_ep] = tree.get_barycentric_data(pts_out,false);

[~,idx_a] = min(D_sp,[],2,'linear');
[~,idx_b] = min(D_ep,[],2,'linear');

% create adjacency matrix from endpoints
node_a = T_sp(idx_a) + int64(hemi_in.' * size(grid.V, 1));
node_b = T_ep(idx_b) + int64(hemi_out.' * size(grid.V, 1));

A = full(sparse(node_a, node_b, 1, 2*size(grid.V,1), 2*size(grid.V,1)));
A = A + A.';

% generate smooth connectome
connectome = kernel * A * kernel.';

end

