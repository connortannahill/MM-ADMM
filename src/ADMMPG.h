#ifndef SOLVER_H
#define SOLVER_H

#include <Eigen/Sparse>
#include <unordered_map>
#include "Assembly.h"
#include <string>

using namespace std;

class ADMMPG {
public:
    ADMMPG(double dt, Assembly &a);
    ~ADMMPG();
    void step(int nIters, double tol);
private:
    Assembly *a;
    Eigen::VectorXd *x;
    Eigen::VectorXd *xPrev;
    Eigen::VectorXd *xBar;
    Eigen::VectorXd *uBar;
    Eigen::VectorXd *z;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> *cgSol;
    Eigen::SparseMatrix<double> *WD_T;
    Eigen::VectorXd *vec;
    Eigen::VectorXd *DXpU;
    // bool matrixFactored;
    double dt;
};

#endif