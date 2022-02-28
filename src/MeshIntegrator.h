#ifndef SOLVER_H
#define SOLVER_H

#include <Eigen/Sparse>
#include <unordered_map>
// #include "Assembly.h"
#include "Mesh.h"
#include <string>

using namespace std;

template <int D>
class MeshIntegrator {
public:
    MeshIntegrator(double dt, Mesh<D> &a);
    ~MeshIntegrator();
    double step(int nIters, double tol);
    double eulerStep(double tol);
    double backwardsEulerStep(double dt, double tol);
    double getEnergy();
    void outputX(const char *fname);
    void outputZ(const char *fname);
    void done();
    double proxTime;
    double multTime;
    double cgTime;
    double predTime;

    Mesh<D> *a;
    Eigen::SparseMatrix<double> *t;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> *cg;
    Eigen::VectorXd *x;
    Eigen::VectorXd *xPrev;
    Eigen::VectorXd *xBar;
    Eigen::VectorXd *uBar;
    Eigen::VectorXd *z;
    Eigen::VectorXd *zPrev;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> *cgSol;
    Eigen::SparseMatrix<double> *WD_T;
    Eigen::SparseMatrix<double> *D_T;
    Eigen::SparseMatrix<double> *WD_TWD;
    Eigen::SparseMatrix<double> *M;
    Eigen::VectorXd *vec;
    Eigen::VectorXd *DXpU;
    // bool matrixFactored;
    double dt;
    double dtPrev;
    double energyCur = INFINITY;
    int stepsTaken;

};

#endif