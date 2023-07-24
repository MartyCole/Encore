function pca(pSublist,pOutput,varargin)

addpath('../interpolation/')
addpath('../tangent_basis')
addpath('../utils/')
addpath('../')

addpath('/nas/longleaf/home/mrcole/SBCI_Toolkit/io')

load(pSublist, 'sub_idx');
load('/nas/longleaf/home/mrcole/ScratchingTheSurface/Code/exampleData/selected_subs_1000_totalcog.mat', 'selected_subs')

% location of the connectivity data we would like to load
input_folder = "/overflow/zzhanglab/SBCI_Finished_ABCD_Data/";
input_subfolder = "/psc_sbci_final_files/sbci_connectome/";
input_file = "mesh_intersections_ico4.mat";

% grid to do evaluations over
ico_grid = load('/nas/longleaf/home/mrcole/ABCD/SBCI_AVG/template_sphere_grid_ico4.mat');
lh_mesh.V = normr(ico_grid.lh_V);
lh_mesh.T = ico_grid.lh_T;
rh_mesh.V = normr(ico_grid.rh_V);
rh_mesh.T = ico_grid.rh_T;

lh_grid = SphericalGrid(lh_mesh,30);
rh_grid = SphericalGrid(rh_mesh,30);

P = length(lh_mesh.V) + length(rh_mesh.V);
N_subs = length(sub_idx);
N_elems = ((P-1)*P) / 2;

elem_idx = logical(triu(ones(P,P), 1));

PCA_mat = zeros(N_subs,N_elems);

for i = 1:N_subs    
    % load the subject endpoint data
    tmp = load(sprintf('%s/%s/%s/%s', input_folder, selected_subs(sub_idx(i)), input_subfolder, input_file));
    
    % get the cartesian coordinates of the endpoints (currently saved as
    % barycentric coordinates, but I will probably change that so we can
    % skip this step in a more finalised version of the pipeline
    [coords_in, coords_out, surf_in, surf_out] = get_endpoint_coords(ico_grid, tmp);
    coords_in = normr(coords_in);
    coords_out = normr(coords_out);

    % initialise the concon object
    concon = Concon(lh_grid,rh_grid,[coords_in,surf_in'],[coords_out,surf_out'],20,0.01,1e-10);
    
    if nargin == 3
        load(sprintf('%s/css_registered_warp%0.4i.mat', varargin{1}, sub_idx(i)),'lh_warp','rh_warp');
        SC = concon.evaluate(lh_warp,rh_warp);
    else
        SC = concon.evaluate();
    end

    PCA_mat(i,:) = log(SC(elem_idx) + 1);
    disp(i)
end

[~,score,latent,explained] = fastpca(PCA_mat);

save(sprintf('%s.mat', pOutput), 'score', 'explained', 'latent', '-v7.3');

end
