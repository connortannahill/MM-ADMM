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
    // Eigen::SparseMatrix<double> t(*(this->a)->M);
    M = new Eigen::SparseMatrix<double>(*(this->a)->M);
    WD_T = new Eigen::SparseMatrix<double>(( *((this->a)->W) * *((this->a)->Dmat)).transpose());
    WD_TWD = new Eigen::SparseMatrix<double>((*WD_T * (*(this->a)->W) * (*(this->a)->Dmat)));

    // *t = *M + dtsq*WD_TWD;
    this->t = new Eigen::SparseMatrix<double>(*M + dtsq*(*WD_TWD));

    // Congugate gradient compute matrix
    this->cg = new Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper>();
    this->cg->compute(*t);

    // Compute the sparse Cholseky factorization.
    cgSol = new Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>>(*t);

    DXpU = new Eigen::VectorXd(*z);
    vec = new Eigen::VectorXd(*z);
}

/**
 * Assumes a constant mass matrix for now
*/
template <int D>
double MeshIntegrator<D>::backwardsEulerStep(double dt, double tol) {
    // cout << "IN BACKWARDS EULER STEP ENERGY INIT = " << a->computeEnergy(*x) << endl;
    Eigen::VectorXd grad(x->size());
    grad.setZero();
    double Ih = a->backwardsEulerStep(dt, *x, grad, tol);

    // cout << "IN BACKWARDS EULER STEP ENERGY FIN = " << a->computeEnergy(*x) << endl;
    return Ih;
}

template <int D>
double MeshIntegrator<D>::getEnergy() {
    return a->computeEnergy(*x);
}

/**
 * Assumes a constant mass matrix for now
*/
template <int D>
double MeshIntegrator<D>::eulerStep(double tol) {
    Eigen::VectorXd grad(x->size());
    double Ih = a->eulerStepMod(*x, grad);
    *x -= (dt / a->tau) * grad;

    // return a->computeEnergy(*x);
    return Ih;
}


/**
 * Assumes a constant mass matrix for now
*/
template <int D>
double MeshIntegrator<D>::step(int nIters, double tol) {

    // Setup the assembly for the step
    (this->a)->setUp();

    // Make prediction for next value of x (sorta) and the next time step
    this->dtPrev = dt;
    double Ih;
    Ih = (this->a)->predictX(dt, energyCur, *this->xPrev, *this->x, *this->xBar);

    // Get xBar, the predicted (explicit) location of the nodes independent of constraints
    double dtsq = dt*dt;

    *xPrev = *x;
    *x = *xBar;
    *z = (*a->Dmat) * *x;
    if (!a->stepTaken) {
        uBar->setZero();
    }

    if (dt != dtPrev) {
        // *t = *M + dtsq*(*WD_TWD);
        // this->cg->compute(*t);
    }

    // Update the solution x^{n+1}
    *vec =  ((a->m) * (*xBar)) + dtsq*(( *WD_T * ((a->w) * (*z - *uBar))));
    // *x = this->cg->solveWithGuess(*vec, *xBar);
    *x = this->cg->solve(*vec);
    // *x = cgSol->solve(*vec);

    int i;
    double IhCur = 0;
    double IhPrev = 0;
    for (i = 0; i < nIters; i++) {
        // Update z_{n+1} using the assembly prox algorithm
        *DXpU = (*(a->Dmat))*(*x) + (*uBar);
        IhCur = a->prox(dt, *x, *DXpU, *z, tol);
        a->stepTaken = true;

        // Update the Lagrange multiplier uBar^{n+1}
        *uBar = *DXpU - *z;

        // Update the solution x^{n+1}
        *vec =  ((a->m) * (*xBar)) + dtsq*(( *WD_T * ((a->w) * (*z - *uBar))));
        // *x = this->cg->solveWithGuess(*vec, *xBar);
        // *x = cgSol->solve(*vec);
        *x = this->cg->solve(*vec);

        if (((*a->Dmat) * *this->x - *z).norm() < tol) {
            break;
        }

        IhPrev = IhCur;
    }
    cout << "ADMM in " << i << " iters" << endl;
    cout << "primal res = " << ((*a->Dmat) * *this->x - *z).norm() << endl;

    // Update the assembly using the new locations
    a->updateAfterStep(dt, *xPrev, *x);

    stepsTaken++;
    // this->energyCur = a->computeEnergy(*x);
    // // return a->eulerStep(*x, *vec);
    // return energyCur;
    return Ih;
}

template <int D>
void MeshIntegrator<D>::done() {
    a->updateAfterStep(dt, *xPrev, *x);
}

template <int D>
MeshIntegrator<D>::~MeshIntegrator() {
    delete x;
    delete xBar;
    delete z;
    delete xPrev;
    // delete cgSol;
    delete WD_T;
    delete DXpU;
    delete vec;

    delete M;
    delete WD_TWD;
    delete cg;

    delete t;
}

template <int D>
void MeshIntegrator<D>::outputX(const char *fname) {
    std::ofstream outFile;
    outFile.open(fname);

    for (int i = 0; i < x->rows()/D; i++) {
        for (int j = 0; j < D-1; j++) {
            outFile << (*x)(i*D+j) << ", ";
        }
        outFile << (*x)(i*D+(D-1)) << endl;
    }

    outFile.close();
}

template <int D>
void MeshIntegrator<D>::outputZ(const char *fname) {
    std::ofstream outFile;
    outFile.open(fname);

    for (int i = 0; i < z->rows()/D; i++) {
        for (int j = 0; j < D-1; j++) {
            outFile << (*z)(i*D+j) << ", ";
        }
        outFile << (*z)(i*D+(D-1)) << endl;
    }

    outFile.close();
}

template class MeshIntegrator<2>;
template class MeshIntegrator<3>;