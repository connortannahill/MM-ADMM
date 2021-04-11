#include <Eigen/Sparse>
#include "MeshIntegrator.h"
#include <vector>
#include <unordered_map>
#include <string>
// #include "Assembly.h"
#include "Mesh.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>

using namespace std;

template <int D>
MeshIntegrator<D>::MeshIntegrator(double dt, Mesh<D> &a) {
    this->a = &a;
    this->dt = dt;

    // Allocate and assign initial values
    x = new Eigen::VectorXd(D*a.getNPnts());
    xPrev = new Eigen::VectorXd(D*a.getNPnts());
    xBar = new Eigen::VectorXd(D*a.getNPnts());

    this->stepsTaken = 0;

    // Assign initial values to x and z (xBar must be assigned at each step)
    a.copyX(*x);
    a.copyX(*xPrev);
    a.copyX(*xBar);

    // Compute the initial z vlaue
    z = new Eigen::VectorXd(*((this->a)->Dmat) * (*this->x));

    uBar = new Eigen::VectorXd(Eigen::VectorXd::Constant(z->size(), 0.0));

    // Prefactor the matrix (assuming constant matrix for now)
    double dtsq = dt*dt;
    Eigen::SparseMatrix<double> t(*(this->a)->M);
    WD_T = new Eigen::SparseMatrix<double>(( *((this->a)->W) * *((this->a)->Dmat)).transpose());

    t = t + dtsq*((*WD_T * (*(this->a)->W) * (*(this->a)->Dmat)));

    // Compute the sparse Cholseky factorization.
    cgSol = new Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>>(t);

    DXpU = new Eigen::VectorXd(*z);
    vec = new Eigen::VectorXd(*z);
}

/**
 * Assumes a constant mass matrix for now
*/
template <int D>
double MeshIntegrator<D>::step(int nIters, double tol) {
    // Get xBar, the predicted (explicit) location of the nodes independent of constraints
    double dtsq = dt*dt;

    // Setup the assembly for the step
    (this->a)->setUp();

    // Make prediction for next value of x (sorta)
    (this->a)->predictX(dt, *this->xPrev, *this->x, *this->xBar);

    *xPrev = *x;
    *z = (*a->Dmat) * *x;

    *x = *xBar;
    // uBar->setZero();

    time_t start, prox, xUpdate, admm;
    prox = 0;
    start = 0;
    xUpdate = 0;
    admm = 0;

    admm = clock();

    int i;
    double IhCur, IhPrev;
    for (i = 0; i < nIters; i++) {
        // cout << "Running prox" << endl;
        // Update z_{n+1} using the assembly prox algorithm
        start = clock();
        *DXpU = (*(a->Dmat))*(*x) + (*uBar);
        // a->updateAfterStep(dt, *xPrev, *x);
        // cout << "running prox" << endl;
        IhCur = a->prox(dt, *x, *DXpU, *z);
        // assert(false);
        // cout << "FINSIEHD running prox" << endl;
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
        // double primalRes = ((*a->Dmat)*(*x) - *z).norm();
        // cout << "Primal res = " << primalRes << endl;
        cout << "check " << abs((IhPrev - IhCur)/(IhPrev)) << endl;
        IhPrev = IhCur;

        string zOutStr = "./gifout/Z"+to_string(stepsTaken)+"-"+to_string(i);
        cout << zOutStr << endl;
        outputZ(zOutStr.c_str());
        if (i >= 1 && abs((IhPrev - IhCur)/(IhPrev)) < tol) {
            break;
        }

    }

    // Update the assembly using the new locations
    start = clock();
    a->updateAfterStep(dt, *xPrev, *x);

    stepsTaken++;
    return IhCur;
}

template <int D>
MeshIntegrator<D>::~MeshIntegrator() {
    delete x;
    delete xBar;
    delete z;
    delete xPrev;
    delete cgSol;
    delete WD_T;
    delete DXpU;
    delete vec;
}

template <int D>
void MeshIntegrator<D>::outputX(const char *fname) {
    std::ofstream outFile;
    outFile.open(fname);
    cout << "outputting x" << endl;
    cout << "size of x (" << x->rows() << ", " << x->cols() << endl;

    for (int i = 0; i < x->rows()/D; i++) {
        for (int j = 0; j < D-1; j++) {
            outFile << (*x)(i*D+j) << ", ";
        }
        outFile << (*x)(i*D+(D-1)) << endl;
    }

    outFile.close();
    cout << "FINSIEHD outputting x" << endl;
}

template <int D>
void MeshIntegrator<D>::outputZ(const char *fname) {
    std::ofstream outFile;
    outFile.open(fname);
    cout << "outputting z" << endl;

    for (int i = 0; i < z->rows()/D; i++) {
        for (int j = 0; j < D-1; j++) {
            outFile << (*z)(i*D+j) << ", ";
        }
        outFile << (*z)(i*D+(D-1)) << endl;
    }

    outFile.close();
    cout << "finished Z" << endl;
}

template class MeshIntegrator<2>;
template class MeshIntegrator<3>;
