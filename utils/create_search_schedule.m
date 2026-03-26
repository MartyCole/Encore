function shells = create_search_schedule()

shells = {
    struct('radius', 20.0, 'dirs', 4, 'K', 5)
    struct('radius', 10.0, 'dirs', 4, 'K', 3)
    struct('radius', 5.0,  'dirs', 5, 'K', 1)
    struct('radius', 2.5,  'dirs', 6, 'K', 1)
    struct('radius', 1.0,  'dirs', 7, 'K', 1)
};

R = cellfun(@(sh) rotation_shell(sh.radius, sh.dirs), shells, 'UniformOutput', false);

for i = 1:5
    shells{i}.R = R{i};
end

end

function S_list = rotation_shell(radius_deg, subdiv)
    % generate directions on an icosphere
    dirs_all = icosphere_dirs(subdiv);

    % select only those within the spherical cap
    dirs = icosphere_cap(dirs_all, radius_deg);

    % preallocate rotation list
    n = size(dirs,1);
    S_list = cell(n+1,1);

    % reference axis (z-axis)
    z_axis = [0;0;1];

    for j = 1:n
        target = dirs(j,:).';

        % compute rotation angle
        cos_theta = dot(z_axis, target);
        cos_theta = max(min(cos_theta,1),-1); 
        theta = acos(cos_theta);

        % if vectors are parallel, rotation is identity
        if theta < 1e-12
            S_list{j} = eye(3);
            continue
        end

        % compute rotation axis
        k = cross(z_axis, target);
        nk = norm(k);

        % safety check
        if nk < 1e-12
            S_list{j} = eye(3);
            continue
        end

        k = k / nk;

        % Rodrigues' rotation formula
        K = [ 0, -k(3), k(2);
              k(3), 0, -k(1);
             -k(2), k(1), 0 ];

        R = eye(3) + sin(theta)*K + (1-cos(theta))*(K*K);

        S_list{j} = R;
    end

    % add identity rotation
    S_list{end} = eye(3);
end

function dirs_cap = icosphere_cap(dirs_all, theta_max_deg)
    theta_max = deg2rad(theta_max_deg);
    z = dirs_all(:,3);
    theta = acos(max(min(z,1),-1));
    dirs_cap = dirs_all(theta <= theta_max, :);
end

function dirs = icosphere_dirs(subdiv)
    ico = icosphere(subdiv);
    dirs = ico.V;
end
