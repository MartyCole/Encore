# Encore - mEtric-based coNtinuous COnnectivity REgistration

Requires libigl for some geometric operation (AABBTree, Voronoi Area, Barycentric Interpolation): https://libigl.github.io/ 

## Overview
TODO: Add description

## Classes
- SphericalGrid - provides geometry, basis, tangent frames
- Kernel / SphericalHeatKernel - provides smoothing kernels + derivatives
- Concon - uses kernel to build continuous connectome
- Encore - registers Concon objects

- SphericalWarp - composed repeatedly during Encore registration

## Pipeline
1. Build spherical grids
2. Construct continuous connectomes
3. Estimate kernel bandwidth
4. Compute template
5. Register subjects
6. Analyze warp fields

## Example

```
ico_mesh = icosphere(ICO_RESOLUTION);  % generate a triagulated mesh for the grid
grid = SphericalGrid(ico_mesh, l);     % generate a grid object from the mesh

% the concon class object representing connectomes
F1 = Concon(grid,grid,F1_start_pts,F1_end_pts,F1_hemi_in,F1_hemi_out);
F2 = Concon(grid,grid,F2_start_pts,F2_end_pts,F2_hemi_in,F2_hemi_out);

% calculate heat kernel for smoothing connectomes
kernel_builder = SphericalHeatKernel(grid,grid,l);

[sigma, ISE] = kernel_builder.cross_validate_sigma(F1,linspace(0.0005,0.005,10));
[K, dK] = kernel_builder.compute_sigma(sigma);

encore = Encore(grid,grid,15,0.05,100,1e-4);
[F1_lh_warp,F1_rh_warp,~] = encore.register(F2,F1,K,dK,'init_rotation',true,'verbose',11);
```

## Citation
TODO: Add reference

