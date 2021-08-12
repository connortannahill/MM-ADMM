#ifndef MESH_H
#define MESH_H

#include "AdaptationFunctional.h"
#include "MonitorFunction.h"
#include "MeshInterpolator.h"
#include "PhaseM.h"
#include <Eigen/Sparse>
#include <vector>
#include <set>

using namespace std;

template <int D=-1>
class Mesh {
public:
    enum NodeType {
        BOUNDARY_FIXED,
        BOUNDARY_FREE,
        INTERIOR
    };
    // Mesh(Eigen::MatrixXd &X, Eigen::MatrixXi &F, Eigen::VectorXi &boundaryMask,
    //         MonitorFunction<D> *M, unordered_map<string, double> params);
    Mesh(Eigen::MatrixXd &X, Eigen::MatrixXi &F, vector<NodeType> &boundaryMask,
            MonitorFunction<D> *M, double rho, double tau);
    void evalMonitorAtPoint(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &mVal);
    void outputSimplices(const char *fname);
    void copyX(Eigen::VectorXd &tar);
    void outputPoints(const char *fname);
    void setUp();
    int getNPnts();
    ~Mesh();
    Eigen::MatrixXd *Vc;
    Eigen::MatrixXd *Vp;
    Eigen::MatrixXi *F;
    Eigen::VectorXd *Ih;
    vector<NodeType> *boundaryMask;
    Eigen::VectorXd *DXpU;
    Eigen::MatrixXd *monitorEvals;
    MonitorFunction<D> *Mon;
    bool stepTaken = false;
    MeshInterpolator<D> *mapEvaluator;
    void outputBoundaryNodes(const char *fname);
    Eigen::FullPivLU<Eigen::Matrix<double, D*(D+1), D*(D+1)>> *lu;

    Eigen::SparseMatrix<double> *M;
    Eigen::SparseMatrix<double> *Dmat;
    Eigen::SparseMatrix<double> *W;

    int a;
    double tau;
    double prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DxpU, Eigen::VectorXd &z);
    void updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x);
    void predictX(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar);
    double newtonOptSimplex(int zId, Eigen::Vector<double, D*(D+1)> &z,
            Eigen::Vector<double, D*(D+1)> &xi, int nIter, double tol);
    void printDiff();
    double BFGSSimplex(int zId, Eigen::Vector<double,D*(D+1)> &z,
        Eigen::Vector<double,D*(D+1)> &xi, int nIter);

    Eigen::MatrixXi *faceList;
    vector<set<int>> *faceConnects;
    void buildFaceList();
    void projectOntoBoundary(int nodeId, Eigen::Vector<double, D> &node);
    void projection2D(int nodeId, Eigen::Vector<double, D> &node);
    void projection3D(int nodeId, Eigen::Vector<double, D> &node);
    
    AdaptationFunctional<D> *I_wx;

    int nx, ny;
    double xa, xb, ya, yb;
    double hx, hy;
    double rho;
    double w;
    int nPnts;

    void buildDMatrix();
    void buildWMatrix(double w);
    void buildMassMatrix(Eigen::VectorXd &m);
};

#endif
