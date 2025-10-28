function [dice,x_dens,y_dens] = get_dice_score(X,Y,lh_grid,rh_grid,threshold)
    lh_aabb = AABBtree(lh_grid);
    rh_aabb = AABBtree(rh_grid);

    lh_P = length(lh_grid.V);
    rh_P = length(rh_grid.V);

    X.st_idx = zeros(size(X.st_points,1),3);
    X.en_idx = zeros(size(X.en_points,1),3);
    Y.st_idx = zeros(size(Y.st_points,1),3);
    Y.en_idx = zeros(size(Y.en_points,1),3);

    X.st_bary = zeros(size(X.st_points,1),3);
    X.en_bary = zeros(size(X.en_points,1),3);
    Y.st_bary = zeros(size(Y.st_points,1),3);
    Y.en_bary = zeros(size(Y.en_points,1),3);

    [X.st_bary(X.st_hemi == 0,:), ...
        X.st_idx(X.st_hemi == 0,:)] = lh_aabb.get_barycentric_data(X.st_points(X.st_hemi == 0,:),false);
    [X.en_bary(X.en_hemi == 0,:), ...
        X.en_idx(X.en_hemi == 0,:)] = lh_aabb.get_barycentric_data(X.en_points(X.en_hemi == 0,:),false);
    [X.st_bary(X.st_hemi == 1,:), ...
        X.st_idx(X.st_hemi == 1,:)] = rh_aabb.get_barycentric_data(X.st_points(X.st_hemi == 1,:),false);
    [X.en_bary(X.en_hemi == 1,:), ...
        X.en_idx(X.en_hemi == 1,:)] = rh_aabb.get_barycentric_data(X.en_points(X.en_hemi == 1,:),false);

    [Y.st_bary(Y.st_hemi == 0,:), ...
        Y.st_idx(Y.st_hemi == 0,:)] = lh_aabb.get_barycentric_data(Y.st_points(Y.st_hemi == 0,:),false);
    [Y.en_bary(Y.en_hemi == 0,:), ...
        Y.en_idx(Y.en_hemi == 0,:)] = lh_aabb.get_barycentric_data(Y.en_points(Y.en_hemi == 0,:),false);
    [Y.st_bary(Y.st_hemi == 1,:), ...
        Y.st_idx(Y.st_hemi == 1,:)] = rh_aabb.get_barycentric_data(Y.st_points(Y.st_hemi == 1,:),false);
    [Y.en_bary(Y.en_hemi == 1,:), ...
        Y.en_idx(Y.en_hemi == 1,:)] = rh_aabb.get_barycentric_data(Y.en_points(Y.en_hemi == 1,:),false);
    
    X.st_idx(X.st_hemi == 1,:) = X.st_idx(X.st_hemi == 1,:) + lh_P;
    X.en_idx(X.en_hemi == 1,:) = X.en_idx(X.en_hemi == 1,:) + lh_P;

    Y.st_idx(Y.st_hemi == 1,:) = Y.st_idx(Y.st_hemi == 1,:) + lh_P;
    Y.en_idx(Y.en_hemi == 1,:) = Y.en_idx(Y.en_hemi == 1,:) + lh_P;

    X.density = accumarray([X.st_idx(:); X.en_idx(:)], [X.st_bary(:); X.en_bary(:)]);
    X.density = X.density ./ (2 * size(X.st_points,1));
    X.density = padarray(X.density,lh_P+rh_P-size(X.density,1),0,'post');

    Y.density = accumarray([Y.st_idx(:); Y.en_idx(:)], [Y.st_bary(:); Y.en_bary(:)]);
    Y.density = Y.density ./ (2 * size(Y.st_points,1));
    Y.density = padarray(Y.density,lh_P+rh_P-size(Y.density,1),0,'post');

    X.set = X.density > threshold;
    Y.set = Y.density > threshold;

    dice = 2*sum(X.set & Y.set) / (sum(X.set) + sum(Y.set));
    x_dens = X.density;
    y_dens = Y.density;
end