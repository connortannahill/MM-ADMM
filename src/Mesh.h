#ifndef MESH_H
#define MESH_H

#include "AdaptationFunctional.h"
#include "MonitorFunction.h"
#include "MeshInterpolator.h"
#include "NodeType.h"
#include <Eigen/Sparse>
#include <vector>
#include <set>

using namespace std;

template <int D=-1>
class Mesh {
public:
    double m;
    // Mesh(Eigen::MatrixXd &X, Eigen::MatrixXi &F, Eigen::VectorXi &boundaryMask,
    //         MonitorFunction<D> *M, unordered_map<string, double> params);
    Mesh(Eigen::MatrixXd &Xc, Eigen::MatrixXd &Xp, Eigen::MatrixXi &F, vector<NodeType> &boundaryMask,
            MonitorFunction<D> *M, int numThreads, double rho, double tau, int integrationMode);
    Mesh(Eigen::MatrixXd &Xp, Eigen::MatrixXi &F, vector<NodeType> &boundaryMask,
            MonitorFunction<D> *M, int numThreads, double rho, double tau, int integratioMode);
    // Mesh(Eigen::MatrixXd &X, Eigen::MatrixXi &F, vector<NodeType> &boundaryMask,
    //         MonitorFunction<D> *M, int numThreads, double rho, double tau);
    void evalMonitorAtPoint(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &mVal);
    void outputSimplices(const char *fname);
    void copyX(Eigen::VectorXd &tar);
    void outputPoints(const char *fname);
    void computeNodalGrads(Eigen::VectorXd &grad);
    void meshInit(Eigen::MatrixXd &Xc, Eigen::MatrixXd &Xp, 
            Eigen::MatrixXi &F, vector<NodeType> &boundaryMask,
            MonitorFunction<D> *Mon, int numThreads, double rho, double tau);
    map<int, vector<int>> *simplexMembership;
    void FSubJac(double dt, int pntId, Eigen::VectorXd &x, Eigen::VectorXd &grad);
//     void FSub(double dt, int sId, Eigen::VectorXd &xk, Eigen::VectorXd &xkp1, Eigen::VectorXd &grad, Eigen::Vector<double, D*(D+1)> &F);
    void buildEulerJac(double dt, Eigen::VectorXd &x, Eigen::VectorXd &grad);
    double backwardsEulerStep(double dt, Eigen::VectorXd &x, Eigen::VectorXd &grad, double tol);
    Eigen::SparseMatrix<double, Eigen::RowMajor> *jac;
    // Eigen::MatrixXd *jac;
    Eigen::SparseMatrix<double, Eigen::RowMajor> *sparseId;
    Eigen::VectorXd *dx;
    Eigen::VectorXd *xn;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > *cg;
    // vector<vector<double*>*> *jacCoefs;
    // Eigen::MatrixXd *jacPrev;
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
    // void sortFaceList();
    MeshInterpolator<D> *mapEvaluator;
    void outputBoundaryNodes(const char *fname);
//     void outputMonitor();
//    Eigen::FullPivLU<Eigen::Matrix<double, D*(D+1), D*(D+1)>> *lu;
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
    double prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DxpU, Eigen::VectorXd &z);
    void reOrientElements(Eigen::MatrixXd &Xp, Eigen::MatrixXi &F);
    void updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x);
    double predictX(double dt, double I, Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar);
    double newtonOptSimplex(int zId, Eigen::Vector<double, D*(D+1)> &z,
            Eigen::Vector<double, D*(D+1)> &xi, int nIter, double tol);
    double bfgsOptSimplex(int zId, Eigen::Vector<double, D*(D+1)> &z,
            Eigen::Vector<double, D*(D+1)> &xi, int nIter, double tol, bool firstStep);
    void printDiff();
    double BFGSSimplex(int zId, Eigen::Vector<double,D*(D+1)> &z,
        Eigen::Vector<double,D*(D+1)> &xi, int nIter);
    double approximateGrads();
    double eulerStep(Eigen::VectorXd &x, Eigen::VectorXd &grad);

    Eigen::MatrixXi *faceList;
    vector<set<int>> *faceConnects;
    vector<int> *pntNeighbours;
    vector<set<int>> *simplexConnects;
    void buildFaceList();
    void buildSimplexMap();
//     void done();
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
