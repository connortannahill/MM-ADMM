#include <Eigen/Sparse>
#include "ADMMPG.h"
#include <vector>
#include <unordered_map>
#include <string>
#include "Assembly.h"
#include <iostream>

using namespace std;

ADMMPG::ADMMPG(double dt, Assembly &a) {
    this->a = &a;
    this->dt = dt;

    // Allocate and assign initial values
    x = new Eigen::VectorXd(a.getD()*a.getNPnts());
    xPrev = new Eigen::VectorXd(a.getD()*a.getNPnts());
    xBar = new Eigen::VectorXd(a.getD()*a.getNPnts());

    // Assign initial values to x and z (xBar must be assigned at each step)
    a.copyX(*x);
    a.copyX(*xPrev);
    a.copyX(*xBar);

    // Compute the initial z vlaue
    z = new Eigen::VectorXd(*((this->a)->D) * (*this->x));

    uBar = new Eigen::VectorXd(Eigen::VectorXd::Constant(z->size(), 0.0));

    // Prefactor the matrix (assuming constant matrix for now)
    double dtsq = dt*dt;
    Eigen::SparseMatrix<double> t(*(this->a)->M);
    WD_T = new Eigen::SparseMatrix<double>(( *((this->a)->W) * *((this->a)->D)).transpose());

    t = t + dtsq*((*WD_T * (*(this->a)->W) * (*(this->a)->D)));

    // Compute the sparse Cholseky factorization.
    cgSol = new Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>>(t);

    DXpU = new Eigen::VectorXd(*z);
    vec = new Eigen::VectorXd(*z);
}

/**
 * Assumes a constant mass matrix for now
*/
void ADMMPG::step(int nIters, double tol) {
    // Get xBar, the predicted (explicit) location of the nodes independent of constraints
    double dtsq = dt*dt;

    // Make prediction for next value of x (sorta)
    (this->a)->predictX(dt, *this->xPrev, *this->x, *this->xBar);

    *xPrev = *x;
    *z = (*a->D) * *x;

    *x = *xBar;
    // uBar->setZero();

    time_t start, prox, xUpdate, admm;
    prox = 0;
    start = 0;
    xUpdate = 0;
    admm = 0;

    admm = clock();

    int i;
    for (i = 0; i < nIters; i++) {
        // cout << "Running prox" << endl;
        // Update z_{n+1} using the assembly prox algorithm
        start = clock();
        *DXpU = (*(a->D))*(*x) + (*uBar);
        a->updateAfterStep(dt, *xPrev, *x);
        a->prox(dt, *x, *DXpU, *z);
        prox += clock() - start;
        // cout << "Finished prox" << endl;

        // Update the Lagrange multiplier uBar^{n+1}
        *uBar = *DXpU - *z;

        // Update the solution x^{n+1}
        start = clock();
        *vec =  ((*a->M) * (*xBar)) + dtsq*(( *WD_T * (*a->W) * (*z - *uBar)));
        *x = cgSol->solve(*vec);
        xUpdate += clock() - start;

        // Compute the primal residual. If it is beneath the tolerance, exit
        double primalRes = ((*a->D)*(*x) - *z).norm();
        if (primalRes < tol) {
            break;
        }
    }

    // cout << "time in full loop " << (clock() - admm)/((double) CLOCKS_PER_SEC) << endl;
    // cout << "time in prox per step " << (prox)/((double) CLOCKS_PER_SEC)/((double)i) << endl;
    // cout << "time in xUpdate per step " << xUpdate/((double) CLOCKS_PER_SEC)/((double)i) << endl;
    // cout << "converged in " << i << " iters" << endl;


    // Update the assembly using the new locations
    start = clock();
    a->updateAfterStep(dt, *xPrev, *x);
    // a->printDiff();
    // cout << "time to update after step " << clock() - start << endl;
}

ADMMPG::~ADMMPG() {
    delete x;
    delete xBar;
    delete z;
    delete xPrev;
    delete cgSol;
    delete WD_T;
    delete DXpU;
    delete vec;
}