#include "mex.h"

#include <Eigen/Core>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/AABB.h>
#include <igl/barycentric_coordinates.h>
#include <igl/slice.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{     
    using namespace Eigen;

    uint64_T* ptr = static_cast<uint64_T*>(mxGetUint64s(prhs[0]));
    igl::AABB<MatrixXd,3>* tree = reinterpret_cast<igl::AABB<MatrixXd,3>*>(*ptr);

    MatrixXd V, Vq, L;    
    MatrixXi F, Fq, I;
    MatrixXd Va, Vb, Vc, C;
    VectorXd sqrD;

    VectorXi xyz(3);
    xyz << 0,1,2;

    igl::matlab::parse_rhs_double(prhs+1,V);
    igl::matlab::parse_rhs_double(prhs+2,F);
    igl::matlab::parse_rhs_double(prhs+3,Vq);  
   
    // find the indices, I, of the intersecting triangles
    tree->squared_distance(V,F,Vq,sqrD,I,C);

    // select the three vertices for each intersecting triangle
    igl::slice(F,I,xyz,Fq);
    igl::slice(V,Fq.col(0),xyz,Va);
    igl::slice(V,Fq.col(1),xyz,Vb);
    igl::slice(V,Fq.col(2),xyz,Vc);
    
    // find the barycentric coordinates, L, for each query point Vq
    igl::barycentric_coordinates(C,Va,Vb,Vc,L);

    // return the barycentric coordinates and intersecting triangles
    igl::matlab::prepare_lhs_double(L,plhs+0);
    igl::matlab::prepare_lhs_double(Fq,plhs+1);
}

/*
mex('query_AABB_tree_mex.cpp', '-R2018a', ...
['-I' fullfile('../../eigen')], ...
['-I' fullfile('../../libigl/include')])
*/