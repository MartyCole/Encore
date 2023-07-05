#include "mex.h"

#include <Eigen/Core>
#include <igl/AABB.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])  
{
    using namespace Eigen;

    uint64_T* ptr = static_cast<uint64_T*>(mxGetUint64s(prhs[0]));
    igl::AABB<MatrixXd,3>* tree = reinterpret_cast<igl::AABB<MatrixXd,3>*>(*ptr);

    delete(tree);
}

/*
mex('destroy_AABB_tree_mex.cpp', '-R2018a', ...
['-I' fullfile('../../eigen')], ...
['-I' fullfile('../../libigl/include')])
*/