function connectome = spherical_kde(grid,kernel,pts_in,pts_out)

% find closet vertex for each point
interpolator = SphericalInterpolator(grid);
[D_sp,T_sp] = interpolator.get_barycentric_data(pts_in);
[D_ep,T_ep] = interpolator.get_barycentric_data(pts_out);

[~,idx_a] = min(D_sp,[],2,'linear');
[~,idx_b] = min(D_ep,[],2,'linear');

% create adjacency matrix from endpoints
node_a = T_sp(idx_a);
node_b = T_ep(idx_b);
A = full(sparse(node_a, node_b, 1, size(grid.V,1), size(grid.V,1)));
A = A + A';

% generate smooth connectome
connectome = kernel * A * kernel';

end

