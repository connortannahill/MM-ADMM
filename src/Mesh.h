#ifndef MESH_H
#define MESH_H

#include "AdaptationFunctional.h"
#include "MonitorFunction.h"
#include "PhaseM.h"
#include <Eigen/Sparse>

template <int D=-1>
class Mesh {
public:
    // Mesh(Eigen::MatrixXd &X, Eigen::MatrixXi &F, Eigen::VectorXi &boundaryMask,
    //         MonitorFunction<D> *M, unordered_map<string, double> params);
    Mesh(Eigen::MatrixXd &X, Eigen::MatrixXi &F, Eigen::VectorXi &boundaryMask,
            PhaseM<D> *M, double rho, double tau);
    void outputSimplices(const char *fname);
    void copyX(Eigen::VectorXd &tar);
    void outputPoints(const char *fname);
    int getNPnts();
    ~Mesh();
    Eigen::MatrixXd *Vc;
    Eigen::MatrixXd *Vp;
    Eigen::MatrixXi *F;
    Eigen::VectorXi *boundaryMask;
    Eigen::VectorXd *DXpU;
    MonitorFunction<D> *Mon;

    Eigen::SparseMatrix<double> *M;
    Eigen::SparseMatrix<double> *Dmat;
    Eigen::SparseMatrix<double> *W;

    int a;
    double tau;
    // void prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DxpU, Eigen::VectorXd &z) override;
    // void updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x) override;
    // void copyX(Eigen::VectorXd &tar) override;
    // void predictX(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar) override;
    void prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DxpU, Eigen::VectorXd &z);
    void updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x);
    void predictX(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar);
    double newtonOptSimplex(int zId, Eigen::Vector<double, D*(D+1)> &z,
            Eigen::Vector<double, D*(D+1)> &xi, int nIter);
    void printDiff();

    AdaptationFunctional<D> *I_wx;

    int nx, ny;
    double xa, xb, ya, yb;
    double hx, hy;
    double rho;
    double w;
    int nPnts;

    // void buildMassMatrix(Eigen::VectorXd &m) override;
    // void buildDMatrix() override;
    // void buildWMatrix(double w) override;
    void buildDMatrix();
    void buildWMatrix(double w);
    void buildMassMatrix(Eigen::VectorXd &m);
};

#endif