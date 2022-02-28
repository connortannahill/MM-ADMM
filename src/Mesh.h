#ifndef MESH_H
#define MESH_H

#include "AdaptationFunctional.h"
#include "MonitorFunction.h"
#include "MeshInterpolator.h"
#include "NodeType.h"
#include "../lib/LASolver/SparseItObj.h"
#include <Eigen/Sparse>
#include <vector>
#include <set>

using namespace std;
using namespace SparseItObj;

template <int D=-1>
class Mesh {
public:
    double m;
    int integrationMode;
    int gradUse;
    Mesh(Eigen::MatrixXd &Xc, Eigen::MatrixXd &Xp, Eigen::MatrixXi &F, vector<NodeType> &boundaryMask,
            MonitorFunction<D> *M, int numThreads, double rho, double w, double tau, int integrationMode, bool gradUse);
    Mesh(Eigen::MatrixXd &Xp, Eigen::MatrixXi &F, vector<NodeType> &boundaryMask,
            MonitorFunction<D> *M, int numThreads, double rho, double w, double tau, int integrationMode, bool gradUse);
    void evalMonitorAtPoint(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &mVal);
    void outputSimplices(const char *fname);
    void copyX(Eigen::VectorXd &tar);
    void outputPoints(const char *fname);
    void computeNodalGrads(Eigen::VectorXd &grad);
    void buildMatrix();
    void meshInit(Eigen::MatrixXd &Xc, Eigen::MatrixXd &Xp, 
            Eigen::MatrixXi &F, vector<NodeType> &boundaryMask,
            MonitorFunction<D> *Mon, int numThreads, double rho, double w, double tau, bool gradUse);
    map<int, vector<int>> *simplexMembership;
    void FSubJac(double dt, int pntId, Eigen::VectorXd &x, Eigen::VectorXd &grad);
    void buildEulerJac(double dt, Eigen::VectorXd &x, Eigen::VectorXd &grad);
    double backwardsEulerStep(double dt, Eigen::VectorXd &x, Eigen::VectorXd &grad, double tol);
    double computeBlockGrad(int zId, Eigen::Vector<double, D*(D+1)> &z,
                Eigen::Vector<double, D*(D+1)> &xi,
                Eigen::Vector<double, D*(D+1)> &grad,
                MeshInterpolator<D> &interp, bool computeGrad,
                bool regularize, double &Igt);

    MatrixIter *jac;
    ParamIter *cgParams;
    double *tol;
    double *rhs;

    Eigen::VectorXd *dx;
    Eigen::VectorXd *xn;

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> *cg;

    void setUp();
    int getNPnts();
    ~Mesh();
    int nSteps = 0;
    Eigen::MatrixXd *Vc;
    Eigen::MatrixXd *Vp;
    Eigen::MatrixXi *F;
    Eigen::VectorXd *Ih;
    vector<NodeType> *boundaryMask;
    Eigen::VectorXd *DXpU;
    Eigen::MatrixXd *monitorEvals;
    MonitorFunction<D> *Mon;
    void setConstant(Eigen::SparseMatrix<double, Eigen::RowMajor> *jac, double a);
    bool stepTaken = false;
    bool hessComputed = false;
    bool compMesh = true;
    MeshInterpolator<D> *mapEvaluator;
    void outputBoundaryNodes(const char *fname);
    vector<Eigen::Matrix<double, D*(D+1), D*(D+1)>> *hessInvs;
    vector<Eigen::Vector<double, D*(D+1)>> *gradCurrs;

    Eigen::SparseMatrix<double> *M;
    Eigen::SparseMatrix<double> *Dmat;
    Eigen::SparseMatrix<double> *W;
    Eigen::VectorXd *grad;
    Eigen::VectorXd *xGlob;
    void fSub (int n, double x[], double fx[], int &iflag );

    int a;
    double tau;
    double prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DxpU, Eigen::VectorXd &z, double tol);
    void reOrientElements(Eigen::MatrixXd &Xp, Eigen::MatrixXi &F);
    void updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x);
    double predictX(double dt, int steps,  Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar);
    double newtonOptSimplex(int zId, Eigen::Vector<double, D*(D+1)> &z,
            Eigen::Vector<double, D*(D+1)> &xi, int nIter, double tol);
    double bfgsOptSimplex(int zId, Eigen::Vector<double, D*(D+1)> &z,
            Eigen::Vector<double, D*(D+1)> &xi, int nIter, double tol, bool firstStep);
    void printDiff();
    double BFGSSimplex(int zId, Eigen::Vector<double,D*(D+1)> &z,
        Eigen::Vector<double,D*(D+1)> &xi, int nIter);
    double approximateGrads();
    double eulerGrad(Eigen::VectorXd &x, Eigen::VectorXd &grad);
    double eulerStep(double dt, Eigen::VectorXd &x, Eigen::VectorXd &grad);
    double eulerStepMod(Eigen::VectorXd &x, Eigen::VectorXd &grad);
    // MatrixIter *matrix;

    Eigen::MatrixXi *faceList;
    vector<set<int>> *faceConnects;
    vector<int> *pntNeighbours;
    vector<set<int>> *simplexConnects;
    void buildFaceList();
    void buildSimplexMap();
    void projectOntoBoundary(int nodeId, Eigen::Vector<double, D> &node);
    void projection2D(int nodeId, Eigen::Vector<double, D> &node);
    void projection3D(int nodeId, Eigen::Vector<double, D> &node);
    double computeEnergy(Eigen::VectorXd &x);
    
    AdaptationFunctional<D> *I_wx;

    int nx, ny;
    double xa, xb, ya, yb;
    double hx, hy;
    double rho;
    double w;
    int nPnts;
    double Ihcur, Ihprev;

    void buildDMatrix();
    void buildWMatrix(double w);
    void buildMassMatrix(Eigen::VectorXd &m);
};

#endif
