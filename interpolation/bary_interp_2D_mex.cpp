#include "mex.h"
#include "matrix.h"

#include <iostream>
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxDouble * lh_points = mxGetDoubles(prhs[0]);
    const mxDouble * rh_points = mxGetDoubles(prhs[1]);
    const mxInt64  * lh_faces  = mxGetInt64s(prhs[2]);
    const mxInt64  * rh_faces  = mxGetInt64s(prhs[3]);

    const mxDouble * data      = mxGetDoubles(prhs[4]);
 
    const int        n_points  = mxGetN(prhs[4]);
    const int        lh_n      = mxGetN(prhs[0]);
    const int        rh_n      = mxGetN(prhs[1]);
    
    const int        n_data    = mxGetNumberOfElements(prhs[4]);

    if (mxGetM(prhs[0]) != 3 && mxGetM(prhs[1]) != 3) {
        std::cout << "Points should be three dimensional" << std::endl;
        return;
    }
    if (mxGetM(prhs[2]) != 3 && mxGetM(prhs[3]) != 3) {
        std::cout << "Faces should be triangular" << std::endl;
        return;
    }
    if (mxGetN(prhs[0]) != mxGetN(prhs[2])  
        && mxGetN(prhs[1]) != mxGetN(prhs[3])) { 
        std::cout << "Inconsistent number of points and faces" << std::endl;
        return;
    }

    const mxDouble * lh_A;
    const mxDouble * rh_A;

    const mxDouble * d0;
    const mxDouble * d1;
    const mxDouble * d2;
 
    const mxInt64  * lh_face;
    const mxInt64  * rh_face;

    double         * new_data;

    size_t out_dims[2];
    out_dims[0] = lh_n;
    out_dims[1] = rh_n;

    plhs[0] = mxCreateNumericArray(2, out_dims, mxDOUBLE_CLASS, mxREAL);    
    new_data = (double *) mxGetDoubles(plhs[0]);
    
    for (int i = 0; i < lh_n; i++) {        
        lh_A    = &lh_points[(i*3)];
        lh_face = &lh_faces[(i*3)];

        d0 = &data[lh_face[0] * n_points];
        d1 = &data[lh_face[1] * n_points];
        d2 = &data[lh_face[2] * n_points];
        
        for (int j = 0; j < rh_n; j++) {
            rh_A    = &rh_points[(j*3)];
            rh_face = &rh_faces[(j*3)];
            
            new_data[(i*lh_n) + j] = (lh_A[0] * d0[rh_face[0]] + lh_A[1] * d1[rh_face[0]] + lh_A[2] * d2[rh_face[0]]) * rh_A[0] +
                                     (lh_A[0] * d0[rh_face[1]] + lh_A[1] * d1[rh_face[1]] + lh_A[2] * d2[rh_face[1]]) * rh_A[1] +
                                     (lh_A[0] * d0[rh_face[2]] + lh_A[1] * d1[rh_face[2]] + lh_A[2] * d2[rh_face[2]]) * rh_A[2];
        }
    }
  
    return;
}

// mex -R2018a bary_interp_2D_mex.cpp