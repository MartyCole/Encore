% FUNCTION NAME:
%   icosphere
%
% DESCRIPTION:
%   Generate a triangle mesh for the unit icosphere in 3 dimensions
%
% INPUT:
%   N - (integer) Recursion level. The mesh will have P = (4^N) * 10 + 2 vertices
%
% OUTPUT:
%   ico - (struct) 
%     ico.V - [Px3] matrix of euclidean coordinates
%     ico.T - [Mx3] matrix of the triangulation 
%
% ASSUMPTIONS AND LIMITATIONS:
%   None
%
function [ico] = icosphere(N)

% sanity check
N = max(min(N,7),0);

% golden ratio
t = (1 + sqrt(5)) / 2;

% twelve vertices of unit icosahedron
ico_vertices = [ -1  +t  +0   % v1
                 +1  +t  +0   % v2
                 -1  -t  +0   % v3
                 +1  -t  +0   % v4
                 +0  -1  +t   % v5
                 +0  +1  +t   % v6
                 +0  -1  -t   % v7
                 +0  +1  -t   % v8
                 +t  +0  -1   % v9
                 +t  +0  +1   % v10
                 -t  +0  -1   % v11
                 -t  +0  +1]; % v12

% normalise vertices to unit size
ico_vertices = normr(ico_vertices);

% faces for unit icosahedron
ico_faces = [ 1  12 6
              1  6  2
              1  2  8
              1  8  11
              1  11 12
              2  6  10
              6  12 5
              12 11 3
              11 8  7
              8  2  9
              4  10 5
              4  5  3
              4  3  7
              4  7  9
              4  9  10
              5  10 6
              3  5  12
              7  3  11
              9  7  8
              10 9  2];
 
% number of final vertices and faces
P = ((4^N) * 10) + 2;
M = (4^N) * 20;

% initialise the mesh
V = zeros(P,3);
T = zeros(M,3);
T(1:20,:) = ico_faces;
V(1:12,:) = ico_vertices;

% index for adding new vertices
n = 13;

% recursively subdivide triangle faces
for level = 1:N           
    Tn = (4^(level-1)) * 20;
    tri = T(1:Tn,:);

    for i = 1:Tn 
        cur_T = tri(i,:);
        
        % calculate mid points of the vertices
        midpoints = sum(cat(3, V(cur_T([1,2,3]),:), V(cur_T([2,3,1]),:)),3).' ./ 2;       

        % find the indices of any new vertices
        match = (V(:,1) == midpoints(1,:)) & ...
                (V(:,2) == midpoints(2,:)) & ...
                (V(:,3) == midpoints(3,:));

        idx = (sum(match) == 0);
        
        % add new points to the triangulation
        n_newpoints = sum(idx);

        if n_newpoints > 0
            V(n:(n+n_newpoints-1),:) = midpoints(:,idx).';
            n = n + n_newpoints;
        end
                
        % store the new triangles              
        idx = (V(:,1) == midpoints(1,:)) & ...
              (V(:,2) == midpoints(2,:)) & ...
              (V(:,3) == midpoints(3,:));

        idx_a = find(idx(:,1));
        idx_b = find(idx(:,2));
        idx_c = find(idx(:,3));

        T(i*4-0,:) = [cur_T(1),idx_a,idx_c];
        T(i*4-1,:) = [idx_a,cur_T(2),idx_b];
        T(i*4-2,:) = [idx_c,idx_b,cur_T(3)];
        T(i*4-3,:) = [idx_a,idx_b,idx_c];             
    end

    % project vertices to the sphere surface
    V = normr(V);
end
    
% return the icosphere mesh
ico = struct('V', V, 'T', T);
end

