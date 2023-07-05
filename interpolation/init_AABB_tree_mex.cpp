#include "mex.h"

#include <Eigen/Core>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/AABB.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{     
    using namespace Eigen;

    MatrixXd V;
    MatrixXi F;

    igl::matlab::parse_rhs_double(prhs+0,V);
    igl::matlab::parse_rhs_double(prhs+1,F);

    igl::AABB<MatrixXd,3>* tree = new igl::AABB<MatrixXd,3>;
    tree->init(V,F);
    
    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    uint64_T* ptr = static_cast<uint64_T*>(mxGetUint64s(plhs[0]));
    *ptr = reinterpret_cast<uint64_T>(tree);   
}

/*
mex('init_AABB_tree_mex.cpp', '-R2018a', ...
['-I' fullfile('../../eigen')], ...
['-I' fullfile('../../libigl/include')])
*/
