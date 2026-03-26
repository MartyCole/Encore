% get_dice_score.m
% Compute the overlap between the density of two sets of endpoints on a mesh grid. 
% Each endpoint is projected onto a triangle in the mesh using barycentric 
% coordinates, and the density per triangle is computed. A Dice score is then 
% calculated based on which triangles exceed a given density threshold.
%
% Syntax:  [dice,x_dens,y_dens] = get_dice_score(X,Y,lh_grid,rh_grid,threshold)
%
% Inputs:
%    X - Endpoints for set A
%    Y - Endpoints for set B
%    lh_grid - Mesh grid for the left hemisphere
%    rh_grid - Mesh grid for the right hemisphere
%    threshold - Only consider triangles with density above threshold
%
% Outputs:
%    dice - The dice score: 2*|X U Y| / (|X| + |Y|);
%    x_dens - The density of X endpoints per triangle
%    y_dens - The density of Y endpoints per triangle
%

function [dice,x_dens,y_dens] = get_dice_score(X,Y,lh_grid,rh_grid,threshold)
    lh_aabb = AABBtree(lh_grid);
    rh_aabb = AABBtree(rh_grid);

    lh_P = length(lh_grid.V);
    rh_P = length(rh_grid.V);

    X = project_and_compute_density(X, lh_aabb, rh_aabb, lh_P, rh_P);
    Y = project_and_compute_density(Y, lh_aabb, rh_aabb, lh_P, rh_P);

    X.set = X.density > threshold;
    Y.set = Y.density > threshold;

    dice = 2*sum(X.set & Y.set) / (sum(X.set) + sum(Y.set));
    x_dens = X.density;
    y_dens = Y.density;
end

function S = project_and_compute_density(S, lh_aabb, rh_aabb, lh_P, rh_P)
    S.st_idx = zeros(size(S.st_points,1),3);
    S.en_idx = zeros(size(S.en_points,1),3);
    S.st_bary = zeros(size(S.st_points,1),3);
    S.en_bary = zeros(size(S.en_points,1),3);

    [S.st_bary(S.st_hemi == 0,:), ...
        S.st_idx(S.st_hemi == 0,:)] = lh_aabb.get_barycentric_data(S.st_points(S.st_hemi == 0,:),false);
    [S.en_bary(S.en_hemi == 0,:), ...
        S.en_idx(S.en_hemi == 0,:)] = lh_aabb.get_barycentric_data(S.en_points(S.en_hemi == 0,:),false);
    [S.st_bary(S.st_hemi == 1,:), ...
        S.st_idx(S.st_hemi == 1,:)] = rh_aabb.get_barycentric_data(S.st_points(S.st_hemi == 1,:),false);
    [S.en_bary(S.en_hemi == 1,:), ...
        S.en_idx(S.en_hemi == 1,:)] = rh_aabb.get_barycentric_data(S.en_points(S.en_hemi == 1,:),false);

    S.st_idx(S.st_hemi == 1,:) = S.st_idx(S.st_hemi == 1,:) + lh_P;
    S.en_idx(S.en_hemi == 1,:) = S.en_idx(S.en_hemi == 1,:) + lh_P;

    S.density = accumarray([S.st_idx(:); S.en_idx(:)], [S.st_bary(:); S.en_bary(:)]);
    S.density = S.density ./ (2 * size(S.st_points,1));
    S.density = padarray(S.density,lh_P+rh_P-size(S.density,1),0,'post');
end
