#include "Assembly.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unordered_map>
#include <iostream>

using namespace std;

template <int D>
Assembly<D>::Assembly() {
    // Allocate the diagonal mass matrix
    // this->nPnts = (int)params["nPnts"];

    dAlloc = false;
    wAlloc = false;
    mAlloc = false;
}

template <int D>
int Assembly<D>::getNPnts() {
    return this->nPnts;
}

template <int D>
void Assembly<D>::setNPnts(int nPnts) {
    this->nPnts = nPnts;
}

template <int D>
void Assembly<D>::buildMassMatrix(Eigen::VectorXd &m) {
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
template <int D>
void Assembly<D>::buildDMatrix() {
    if (!dAlloc) {
        Dmat = new Eigen::SparseMatrix<double>(D*nPnts, D*nPnts);
        dAlloc = true;
    } else {
        assert(false);
    }

    Dmat->reserve(Eigen::VectorXd::Constant(D*nPnts, 1));

    for (int i = 0; i < Dmat->rows(); i++) {
        Dmat->insert(i, i) = 1.0;
    }
}

/**
 * Default implementation, W is the identity matrix
*/
template <int D>
void Assembly<D>::buildWMatrix(double w) {
    this->w = w;

    if (!wAlloc) {
        W = new Eigen::SparseMatrix<double>(D*nPnts, D*nPnts);
        wAlloc = true;
    } else {
        assert(false);
    }

    W->reserve(Eigen::VectorXd::Constant(D*nPnts, 1));

    for (int i = 0; i < W->rows(); i++) {
        W->insert(i, i) = w;
    }
}

template <int D>
Assembly<D>::~Assembly() {
    if (mAlloc) {
        delete M;
    }
    if (dAlloc) {
        delete Dmat;
    }
    if (wAlloc) {
        delete W;
    }
}

template class Assembly<2>;
template class Assembly<3>;