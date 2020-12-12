#include "Assembly.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unordered_map>
#include <iostream>

using namespace std;

Assembly::Assembly(unordered_map<string, double> params) {
    // Allocate the diagonal mass matrix
    this->nPnts = (int)params["nPnts"];

    this->d = (int)params["d"];
    M = new Eigen::SparseMatrix<double>(d*nPnts, d*nPnts);

    dAlloc = false;
    wAlloc = false;
    mAlloc = false;
}

// void Assembly::printDiff() {
//     cout << "Difference after step: " << (*Vp - *Vc).norm() << endl;
// }

int Assembly::getD() {
    return this->d;
}

int Assembly::getNPnts() {
    return this->nPnts;
}
void Assembly::setNPnts(int nPnts) {
    this->nPnts = nPnts;
}
void Assembly::setD(int d) {
    this->d = d;
}

void Assembly::buildMassMatrix(Eigen::VectorXd &m) {
    if (!mAlloc) {
        M = new Eigen::SparseMatrix<double>(m.size(), m.size());
        mAlloc = true;
    } else {
        assert(false);
    }
    M->reserve(Eigen::VectorXd::Constant(m.size(), 1));

    for (int i = 0; i < m.size(); i++) {
        M->insert(i, i) = m[i];
    }
}

/**
 * Default implementation, D is the identity matrix
*/
void Assembly::buildDMatrix() {
    if (!dAlloc) {
        D = new Eigen::SparseMatrix<double>(d*nPnts, d*nPnts);
        dAlloc = true;
    } else {
        assert(false);
    }

    D->reserve(Eigen::VectorXd::Constant(d*nPnts, 1));

    for (int i = 0; i < D->rows(); i++) {
        D->insert(i, i) = 1.0;
    }
}

/**
 * Default implementation, W is the identity matrix
*/
void Assembly::buildWMatrix(double w) {
    this->w = w;

    if (!wAlloc) {
        W = new Eigen::SparseMatrix<double>(d*nPnts, d*nPnts);
        wAlloc = true;
    } else {
        assert(false);
    }

    W->reserve(Eigen::VectorXd::Constant(d*nPnts, 1));

    for (int i = 0; i < W->rows(); i++) {
        W->insert(i, i) = w;
    }
}

Assembly::~Assembly() {
    delete M;
    if (dAlloc) {
        delete D;
    }
    if (wAlloc) {
        delete W;
    }
}