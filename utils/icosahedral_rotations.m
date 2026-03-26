function Rset = icosahedral_rotations()
    % return the 60 rotation matrices of the icosahedral group (SO(3))
    
    % rotation generators
    f = (1 + sqrt(5)) / 2;

    R3 = (1/2) * [1-f,f,-1;-f,-1,1-f;-1,f-1,f];  
    R2 = [-1,0,0;0,-1,0;0,0,1];

    matrices = {R3,R2}; 
   
    idx=1;

    % generate N-dimensional grid of indices
    N = 12;
    v = [1,2];
    grids = cell(1, N);
    [grids{:}] = ndgrid(1:numel(v));

    % flatten and map indices to values in v
    combos = zeros(numel(grids{1}), N);
    for k = 1:N
        combos(:,k) = v(grids{k}(:));
    end

    Rset = cell(size(combos, 1),1);

    for j = 1:size(combos, 1)
        result = eye(3);
        for k = 1:N
            result = round(result * matrices{combos(j,k)},6);
        end
    
        Rset{idx} = result;
        idx = idx + 1;
    end        

    [~,idx] = unique(cellfun(@(x) mat2str(round(x,6)), Rset, 'UniformOutput', false));
    
    Rset = reshape(cell2mat(Rset(idx))',3,3,[]);
end
