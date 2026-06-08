# Example Analysis Pipeline  
### Using 10 Participants from the HCP Young Adult Dataset

The file `example_analysis.m` implements the workflow described in this document.  
The script performs a full **population‑level registration** of structural connectomes using spherical representations of streamline endpoints.

The pipeline consists of:

1. **Build a spherical grid**  
2. **Load subject‑level streamline endpoint data**  
3. **Construct structural connectome objects (`SConcon`)**  
4. **Build a smoothing kernel on the sphere**  
5. **Estimate an optimal diffusion parameter (`sigma`)**  
6. **Compute the kernel and its derivative**  
7. **Estimate a population template**  
8. **Register each subject to the template**  
9. **Save the resulting warp fields**

The output is a set of **diffeomorphic warp fields** for each hemisphere of each subject.

## 1. Grid Construction

```matlab
ICO_RESOLUTION = 4;
l = 30;

ico_mesh = icosphere(ICO_RESOLUTION);
grid = SphericalGrid(ico_mesh, l);
```

### Description

- Generates an **icosphere mesh** at the specified resolution  
- Constructs a **SphericalGrid** object that:
  - Defines sampling points on the sphere  
  - Precomputes spherical harmonic bases up to order `l`  
  - Supports interpolation, smoothing, and registration  

In this example, the same grid is used for both hemispheres, though the framework supports using different grids if needed.

## 2. Loading Subject Data

```matlab
sublist = readlines("../data/sublist", "EmptyLineRule", "skip");
N = length(sublist);
Fs = cell(N,1);
```

Each subject has a corresponding file:

```
<subject>_reg_sphere_intersections.mat
```

containing:

- `surf_in`, `surf_out` - hemisphere labels  
- `vtx_in`, `vtx_out` - streamline endpoints on the sphere  

These are used to construct a **structural continuous connectivity object**:

```matlab
Fs{i} = SConcon(grid,grid,F1_start_pts,F1_end_pts,F1_hemi_in,F1_hemi_out);
```

### Description

`SConcon` represents structural connectivity using streamline endpoints:

- Start point on the sphere  
- End point on the sphere  
- Hemisphere membership  
- Mapping to the spherical grid  

The object supports smoothing, kernel evaluation, and registration.

## 3. Building the Smoothing Kernel

```matlab
kernel_builder = SphericalHeatKernel(grid,grid,l);
[sigma, LOO] = kernel_builder.cross_validate_sigma(Fs{2}, linspace(0.001,0.008,22));
[K, dK] = kernel_builder.compute_sigma(0.005);
```

### Description

The **SphericalHeatKernel** applies heat‑diffusion smoothing to spherical connectomes.

- `cross_validate_sigma` estimates an optimal diffusion time  
- `compute_sigma` evaluates:
  - `K` - the smoothing kernel  
  - `dK` - its derivative (required for registration gradients)

Although each subject may have a different optimal bandwidth, the analysis uses a fixed  
`sigma = 0.005` to ensure consistent smoothing across subjects.

## 4. Registration and Template Estimation

```matlab
encore = Encore(grid,grid,15,0.05,100,1e-4);
template = encore.get_template(Fs,K,100);
```

### The `Encore` Object

`Encore` implements a **diffeomorphic registration algorithm** on the sphere.

Parameters:

- `15` - spherical harmonic order for the velocity field  
- `0.05` - gradient‑descent step size  
- `100` - maximum iterations  
- `1e‑4` - convergence tolerance  

### Template Estimation

`get_template` computes a **population‑average connectome** by iteratively:

1. Computing the median connectome (evaluated using kernel `K`)  
2. Updating the estimate toward the Karcher mean  
3. Refining the template  
4. Repeating until convergence  

## 5. Registering Each Subject

```matlab
[lh_warp, rh_warp, ~] = encore.register(template, Fs{i}, K, dK, "verbose", 1);
save(sprintf("../example_results/%s_warp.mat",sublist(i)), "lh_warp", "rh_warp")
```

For each subject:

- Computes a left‑hemisphere warp  
- Computes a right‑hemisphere warp  
- Saves both to `example_results/`

### Output Format

Each warp contains:

- A diffeomorphic mapping from subject space (originally fsaverage) to template space  
- Represented as **stationary velocity fields (SVFs)**  
- One warp per hemisphere  

These warps can be used to:

- Align connectomes  
- Compare structural connectivity across subjects  
- Visualize deformation fields  
- Perform group‑level statistical analyses  

## 6. How to Run the Example Pipeline

### Requirements

- MATLAB R2022b or later  
- Project codebase located at `../../` relative to the script  
- Data stored in `../data/`  
- A `sublist` file listing subject IDs  

### Running the Script

```matlab
cd path/to/script
example_analysis
```

Results will be written to:

```
../example_results/
```

## 7. Modifying the Pipeline

### Change Grid Resolution

```matlab
ICO_RESOLUTION = 5;
```

Higher resolution -> more vertices -> slower but potentially more precise.

### Change Smoothness of Diffeomorphisms

```matlab
encore = Encore(grid,grid,20,0.05,100,1e-4);
```

Lower spherical harmonic order -> smoother but less expressive warps.

### Use a Different Sigma

```matlab
[K, dK] = kernel_builder.compute_sigma(sigma);
```

# Analyzing Warped Connectomes

After generating warp fields, they can be applied to the `SConcon` objects for analysis.

## Applying a Warp

```matlab
load(sprintf("../example_results/%s_warp.mat",sublist(i)), "lh_warp", "rh_warp")
Fs{i}.warp_connectome(lh_warp, rh_warp);
concon = gather(Fs{i}.evaluate(K));
```

`gather` is required because enCoRe uses the GPU by default.

The resulting matrix `concon` is a **P × P continuous structural connectivity matrix**.  
For numerical stability, it is often helpful to scale the matrix (e.g., multiply by `1e6`).

The repository includes:

- `get_dice_score`  
- `get_overlap_coefficient`  

which can be used to reproduce the evaluation metrics from the paper.