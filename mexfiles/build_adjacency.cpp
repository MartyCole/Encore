#include "mex.h"
#include <Eigen/Dense>

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Vector3f;
using Eigen::Matrix3f;

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6)
        mexErrMsgTxt("Usage: adjacency = build_adjacency_fast(st_idx, en_idx, st_points, en_points, nRows, nCols)");

    const int32_T *st_idx = (const int32_T*) mxGetData(prhs[0]);
    const int32_T *en_idx = (const int32_T*) mxGetData(prhs[1]);
    const double  *st_pts = mxGetPr(prhs[2]);
    const double  *en_pts = mxGetPr(prhs[3]);

    mwIndex N     = mxGetM(prhs[0]);
    mwIndex nRows = (mwIndex) mxGetScalar(prhs[4]);
    mwIndex nCols = (mwIndex) mxGetScalar(prhs[5]);

    plhs[0] = mxCreateNumericMatrix(nRows, nCols, mxSINGLE_CLASS, mxREAL);
    float *adj_ptr = (float*) mxGetData(plhs[0]);

    Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> adj(adj_ptr, nRows, nCols);
    adj.setZero();

    mwIndex ld = N;

    const int32_T *st0 = st_idx;
    const int32_T *st1 = st_idx + ld;
    const int32_T *st2 = st_idx + 2 * ld;

    const int32_T *en0 = en_idx;
    const int32_T *en1 = en_idx + ld;
    const int32_T *en2 = en_idx + 2 * ld;

    const double *sp0 = st_pts;
    const double *sp1 = st_pts + ld;
    const double *sp2 = st_pts + 2 * ld;

    const double *ep0 = en_pts;
    const double *ep1 = en_pts + ld;
    const double *ep2 = en_pts + 2 * ld;

    for (mwIndex k = 0; k < N; k++) {
        Vector3f s(sp0[k], sp1[k], sp2[k]);
        Vector3f e(ep0[k], ep1[k], ep2[k]);

        Matrix3f block = s * e.transpose();   // vectorized outer product

        int r0 = st0[k] - 1;
        int r1 = st1[k] - 1;
        int r2 = st2[k] - 1;

        int c0 = en0[k] - 1;
        int c1 = en1[k] - 1;
        int c2 = en2[k] - 1;

        // scatter-add the 3×3 block
        adj(r0, c0) += block(0,0);
        adj(r0, c1) += block(0,1);
        adj(r0, c2) += block(0,2);

        adj(r1, c0) += block(1,0);
        adj(r1, c1) += block(1,1);
        adj(r1, c2) += block(1,2);

        adj(r2, c0) += block(2,0);
        adj(r2, c1) += block(2,1);
        adj(r2, c2) += block(2,2);
    }
}

// mex -R2018a CXXFLAGS="\$CXXFLAGS -O3 -march=native" -I../../eigen build_adjacency.cpp

