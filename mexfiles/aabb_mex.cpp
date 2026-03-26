#include <Eigen/Core>
#include <igl/AABB.h>
#include <igl/barycentric_coordinates.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>

#include <map>
#include <cstring>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;

struct TreeData {
    igl::AABB<MatrixXd,3>* tree = nullptr;
    MatrixXd V;
    MatrixXi F;

    // reusable buffers for queries
    // resizing is faster than allocating
    VectorXd sqrD;
    VectorXi I;
    MatrixXd C;
    MatrixXi Fq;
    MatrixXd Va, Vb, Vc;
    MatrixXd L;
};

static std::map<int, TreeData> g_trees;
static int g_next_id = 1;

void cleanup()
{
    for (auto& kv : g_trees) {
        delete kv.second.tree;
        kv.second.tree = nullptr;
    }
    g_trees.clear();
}

void cmd_init(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
    if (nrhs != 3)
        mexErrMsgTxt("Usage: id = aabb_mex('init', V, F)");

    TreeData data;

    igl::matlab::parse_rhs_double(prhs+1, data.V);
    igl::matlab::parse_rhs_double(prhs+2, data.F);

    data.tree = new igl::AABB<MatrixXd,3>;
    data.tree->init(data.V, data.F);

    int id = g_next_id++;
    g_trees[id] = std::move(data);

    plhs[0] = mxCreateDoubleScalar(static_cast<double>(id));
}

void cmd_query(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
    if (nrhs != 3)
        mexErrMsgTxt("Usage: [L, Fq] = aabb_mex('query', id, Vq)");

    int id = static_cast<int>(mxGetScalar(prhs[1]));
    auto it = g_trees.find(id);
    if (it == g_trees.end())
        mexErrMsgTxt("Invalid tree ID");

    TreeData& data = it->second;

    MatrixXd Vq;
    igl::matlab::parse_rhs_double(prhs+2, Vq);

    const int m = static_cast<int>(Vq.rows());

    // Resize reusable buffers
    data.sqrD.resize(m);
    data.I.resize(m);
    data.C.resize(m, 3);

    // Query AABB tree
    data.tree->squared_distance(data.V, data.F, Vq,
                                data.sqrD, data.I, data.C);

    // Gather faces (Fq) using indexed access
    data.Fq.resize(m, 3);
    data.Fq = data.F(data.I, Eigen::all);

    // Gather triangle vertices Va, Vb, Vc
    data.Va.resize(m, 3);
    data.Vb.resize(m, 3);
    data.Vc.resize(m, 3);

    data.Va = data.V(data.Fq.col(0), Eigen::all);
    data.Vb = data.V(data.Fq.col(1), Eigen::all);
    data.Vc = data.V(data.Fq.col(2), Eigen::all);

    // Compute barycentric coordinates
    igl::barycentric_coordinates(data.C, data.Va, data.Vb, data.Vc, data.L);
    
    // Return L and Fq
    igl::matlab::prepare_lhs_double(data.L,  plhs+0);
    igl::matlab::prepare_lhs_double(data.Fq, plhs+1);
}

void cmd_destroy(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
    if (nrhs != 2)
        mexErrMsgTxt("Usage: aabb_mex('destroy', id)");

    int id = static_cast<int>(mxGetScalar(prhs[1]));
    auto it = g_trees.find(id);
    if (it == g_trees.end())
        mexErrMsgTxt("Invalid tree ID");

    delete it->second.tree;
    it->second.tree = nullptr;
    g_trees.erase(it);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1 || !mxIsChar(prhs[0]))
        mexErrMsgTxt("First argument must be a command string");

    static bool registered = false;
    if (!registered) {
        mexAtExit(cleanup);
        registered = true;
    }

    char cmd[32];
    mxGetString(prhs[0], cmd, sizeof(cmd));

    if (!std::strcmp(cmd, "init")) {
        cmd_init(nlhs, plhs, nrhs, prhs);
    } else if (!std::strcmp(cmd, "query")) {
        cmd_query(nlhs, plhs, nrhs, prhs);
    } else if (!std::strcmp(cmd, "destroy")) {
        cmd_destroy(nlhs, plhs, nrhs, prhs);
    } else {
        mexErrMsgTxt("Unknown command");
    }
}

// mex -R2018a CXX=g++-9 CC=gcc-9 CXXFLAGS="\$CXXFLAGS -O3 -march=native -ffast-math" -I../../eigen -I../../libigl/include aabb_mex.cpp
