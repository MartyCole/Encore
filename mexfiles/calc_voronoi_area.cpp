#include "mex.h"

#include <Eigen/Core>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/massmatrix.h>
#include <igl/diag.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  using namespace Eigen;

  if (nrhs != 2) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
        "calc_voronoi_area requires 2 input arguments, vertices and faces");
  }
 
  MatrixXd V, D;
  MatrixXi F;
  SparseMatrix<double> M;
  
  igl::matlab::parse_rhs_double(prhs+0,V);
  igl::matlab::parse_rhs_double(prhs+1,F);
   
  F.array() -= 1;

  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
  igl::diag(M,D);

  igl::matlab::prepare_lhs_double(D,plhs); 

  return;
}
/*
mex('calc_voronoi_area.cpp', ... ...
    ['-I' fullfile('../eigen')], ...
    ['-I' fullfile('../libigl/include')])
*/
