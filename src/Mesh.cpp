#include "Mesh.h"
#include "Assembly.h"
#include <Eigen/Dense>
#include <unordered_map>
#include "MonitorFunction.h"
#include "HuangFunctional.h"
#include <stdlib.h>
#include <string>
#include <map>
#include <iostream>
// #include <LBFGS.h>
#include <fstream>
#include <omp.h>
#include <unistd.h>
#include <algorithm>

#define THREAD_NUM 1

using namespace std;
// using namespace LBFGSpp;


void segmentDS(Eigen::VectorXd &x, int xOff, Eigen::Vector2d &z, int num) {
    for (int i = 0; i < num; i++) {
        z(i) = x(xOff+i);
    }
}

void segmentDD(Eigen::VectorXd &x, int xOff, Eigen::VectorXd &z, int zOff, int num) {
    for (int i = 0; i < num; i++) {
        z(zOff+i) = x(xOff+i);
    }
}

void segmentSD(Eigen::Vector2d &z, Eigen::VectorXd &x, int xOff, int num) {
    for (int i = 0; i < num; i++) {
        x(xOff+i) = z(i);
    }
}

template <int D>
int Mesh<D>::getNPnts() {
    return this->nPnts;
}

template <int D>
Mesh<D>::Mesh(Eigen::MatrixXd &X, Eigen::MatrixXi &F, Eigen::VectorXi &boundaryMask,
            MonitorFunction<D> *Mon, double rho, double tau) {

    this->Vc = &X;
    this->Vp = new Eigen::MatrixXd(*this->Vc);
    this->F = &F;
    this->boundaryMask = &boundaryMask;

    this->Ih = new Eigen::VectorXd(F.rows());

    //omp_set_thread_num(THREAD_NUM);
    omp_set_num_threads(THREAD_NUM);

    // Create mesh interpolator
    cout << "Creating the mesh interpolator" << endl;
    mapEvaluator = new MeshInterpolator<D>();
    mapEvaluator->updateMesh((*this->Vc), (*this->F));
    mapEvaluator->interpolateMonitor(*Mon);
    cout << "FINSIEHD Creating the mesh interpolator" << endl;

    this->nPnts = X.rows();

    this->Mon = Mon;

    // Trivial edge matrix
    int nPnts = Vc->rows();

    // Build the monitor simplex values
    monitorEvals = new Eigen::MatrixXd(X.rows(), D*D);

    // Build the mass matrix
    cout << "bulding the mass matrix" << endl;
    Eigen::VectorXd m(Eigen::VectorXd::Constant(nPnts*D, tau));
    buildMassMatrix(m);
    cout << "FINISHED bulding the mass matrix" << endl;

    cout << "bulding the D matrix" << endl;
    buildDMatrix();
    this->w = sqrt(rho);
    buildWMatrix(this->w);
    cout << "FINISHED bulding the D matrix" << endl;

    DXpU = new Eigen::VectorXd((*Dmat)*(m));

    // Create the functional for each vertex
    cout << "bulding the functional" << endl;
    I_wx = new HuangFunctional<D>(*Vc, *Vp, *(this->F), *DXpU, Mon, this->w, 0.0, 0.0);
    cout << "FINSIEHD bulding the functional" << endl;
}

// TODO: add external forces (would be handled in this object, as the meshing object should control its own props)
template <int D>
void Mesh<D>::predictX(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar) {
    xBar = 2.0*x - xPrev;
}

template <int D>
void Mesh<D>::buildMassMatrix(Eigen::VectorXd &m) {
    M = new Eigen::SparseMatrix<double>(m.size(), m.size());

    M->reserve(Eigen::VectorXd::Constant(m.size(), 1));

    for (int i = 0; i < m.size(); i++) {
        M->insert(i, i) = m[i];
    }
}

/**
 * Default implementation, W is the identity matrix
*/
template <int D>
void Mesh<D>::buildWMatrix(double w) {
    this->w = w;

    W = new Eigen::SparseMatrix<double>(F->rows()*D*(D+1), F->rows()*D*(D+1));

    W->reserve(Eigen::VectorXd::Constant(F->rows()*D*(D+1), 1));

    for (int i = 0; i < W->rows(); i++) {
        W->insert(i, i) = w;
    }
}

template <int D>
void Mesh<D>::buildDMatrix() {
    // // Temp matrix (we insert the transpose) and we keep it this way for storage reasons
    Eigen::SparseMatrix<double> Dt(D*Vp->rows(), D*(D+1)*F->rows());

    // // Reserve the memory
    // Dt.reserve(Eigen::VectorXi::Constant(D*(D+1)*F->rows(), 1));

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(D*(D+1)*F->rows());
    // for(...)
    // {
    // // ...
    // tripletList.push_back(T(i,j,v_ij));
    // }
    // SparseMatrixType mat(rows,cols);
    // mat.setFromTriplets(tripletList.begin(), tripletList.end());
    // mat is ready to go!


    // Add the spring constraints to the system
    int pntIds[3];
    int rowStart0, rowStart1, rowStart2, colStart0, colStart1, colStart2;
    // Note: the is columns in the matrix transposed
    for (int col = 0; col < F->rows(); col++) {
        // Extract the ids
        pntIds[0] = (*F)(col, 0);
        pntIds[1] = (*F)(col, 1);
        pntIds[2] = (*F)(col, 2);

        colStart0 = D*(D+1)*col;
        colStart1 = colStart0+D;
        colStart2 = colStart1+D;

        rowStart0 = pntIds[0]*D;
        rowStart1 = pntIds[1]*D;
        rowStart2 = pntIds[2]*D;

        // First block
        for (int n = 0; n < D; n++) {
            for (int m = 0; m < D; m++) {
                if (n == m) {
                    tripletList.push_back(T(n+rowStart0, m+colStart0, 1.0));
                    tripletList.push_back(T(n+rowStart1, m+colStart1, 1.0));
                    tripletList.push_back(T(n+rowStart2, m+colStart2, 1.0));

                    // Dt.insert(n+rowStart0, m+colStart0) = 1.0;
                    // Dt.insert(n+rowStart1, m+colStart1) = 1.0;
                    // Dt.insert(n+rowStart2, m+colStart2) = 1.0;
                } else {
                    tripletList.push_back(T(n+rowStart0, m+colStart0, 0.0));
                    tripletList.push_back(T(n+rowStart1, m+colStart1, 0.0));
                    tripletList.push_back(T(n+rowStart2, m+colStart2, 0.0));
                    // Dt.insert(n+rowStart0, m+colStart0) = 0.0;
                    // Dt.insert(n+rowStart1, m+colStart1) = 0.0;
                    // Dt.insert(n+rowStart2, m+colStart2) = 0.0;
                }
            }
        }
    }

    Dt.setFromTriplets(tripletList.begin(), tripletList.end());

    // Create the D matrix
    Dmat = new Eigen::SparseMatrix<double>(Dt.transpose());

    /**
    // Temp matrix (we insert the transpose) and we keep it this way for storage reasons
    Eigen::SparseMatrix<double> Dt(D*Vp->rows(), D*(D+1)*F->rows());

    // Reserve the memory
    cout << "reserving the memory" << endl;
    Dt.reserve(Eigen::VectorXi::Constant(D*(D+1)*F->rows(), 1));
    cout << "FINISEHD reserving the memory" << endl;

    // Add the spring constraints to the system
    int pntIds[3];
    int rowStart0, rowStart1, rowStart2, colStart0, colStart1, colStart2;
    // Note: the is columns in the matrix transposed
    for (int col = 0; col < F->rows(); col++) {
	cout << "col " << col << " out of " << F->rows() << endl;
        // Extract the ids
        pntIds[0] = (*F)(col, 0);
        pntIds[1] = (*F)(col, 1);
        pntIds[2] = (*F)(col, 2);

        colStart0 = D*(D+1)*col;
        colStart1 = colStart0+D;
        colStart2 = colStart1+D;

        rowStart0 = pntIds[0]*D;
        rowStart1 = pntIds[1]*D;
        rowStart2 = pntIds[2]*D;

        // First block
        for (int n = 0; n < D; n++) {
            for (int m = 0; m < D; m++) {
                if (n == m) {
                    Dt.insert(n+rowStart0, m+colStart0) = 1.0;
                    Dt.insert(n+rowStart1, m+colStart1) = 1.0;
                    Dt.insert(n+rowStart2, m+colStart2) = 1.0;
                } else {
                    Dt.insert(n+rowStart0, m+colStart0) = 0.0;
                    Dt.insert(n+rowStart1, m+colStart1) = 0.0;
                    Dt.insert(n+rowStart2, m+colStart2) = 0.0;
                }
            }
        }
    }

    // Create the D matrix

    Dmat = new Eigen::SparseMatrix<double>(Dt.transpose());
    */

}

/**
 * BFGS over a input simplex
*/
template <int D>
double Mesh<D>::BFGSSimplex(int zId, Eigen::Vector<double,D*(D+1)> &z,
        Eigen::Vector<double,D*(D+1)> &xi, int nIter) {
    return 0;
}

/**
 * Newton's method over a single simplex
*/
template <int D>
double Mesh<D>::newtonOptSimplex(int zId, Eigen::Vector<double, D*(D+1)> &z,
        Eigen::Vector<double, D*(D+1)> &xi, int nIter) {
    double h = 2.0*sqrt(std::numeric_limits<double>::epsilon());
    // double h = 1e-10;
    const int MAX_LS = 10;

    Eigen::Vector<double, D*(D+1)> zPurt;
    Eigen::Vector<double, D*(D+1)> gradZ;
    Eigen::Vector<double, D*(D+1)> gradZPurt;
    Eigen::Matrix<double, D*(D+1), D*(D+1)> hess;
    Eigen::Vector<double, D*(D+1)> p;
    Eigen::Vector<double, D*(D+1)> gradTemp;

    double Ix;
    double Ipurt;

    for (int iters = 0; iters < nIter; iters++) {
        // cout << "Running block grad" << endl;
        Ix = I_wx->blockGrad(zId, z, xi, gradZ, *mapEvaluator);
        // assert(false);

        zPurt = z;

        // Compute the Hessian column-wise
        for (int i = 0; i < D*(D+1); i++) {
            // Compute purturbation
            zPurt(i) += h;

            // Compute gradient at purturbed point
            Ipurt = I_wx->blockGrad(zId, zPurt, xi, gradZPurt, *mapEvaluator);

            hess.col(i) = (gradZPurt - gradZ)/h;

            zPurt(i) = z(i);
        }

        // Compute the Newton direction
        // p = hess.colPivHouseholderQr().solve(-gradZ);
        p = hess.lu().solve(-gradZ);
        zPurt = z + p;

        // Perform backtracking line search in the Hessian direction (should work if Hess is pos def)
        // int lsIters = 0;
        // double alpha = 1.0;
        // double c1 = 0.0;
        // double c2 = 0.9;
        // Ipurt = I_wx->blockGrad(zId, zPurt, xi, gradTemp, *mapEvaluator);

        // while (Ipurt >= Ix - c1*alpha*(gradZ.dot(p)) &&
        //         -p.dot(gradTemp) >= -c2*p.dot(gradZ)  && lsIters < MAX_LS) {
        //     alpha /= 10.0;

        //     zPurt = z + alpha*p;
        //     Ipurt = I_wx->blockGrad(zId, zPurt, xi, gradTemp, *mapEvaluator);

        //     lsIters++;
        // }

        // if (lsIters == MAX_LS) {
        //     break;
        // } else {
            z = zPurt;
        // }
    }
    // assert(false);

    return Ix;

}

template <int D>
double Mesh<D>::prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DXpU, Eigen::VectorXd &z) {
    // Copy DXpU address to local pointer
    *this->DXpU = DXpU;

    mapEvaluator->interpolateMonitor(*this->Mon);
    // mapEvaluator->outputStuff();
    // assert(false);

    // Run Newton's method on each simplex
  //  Eigen::Vector<double, D*(D+1)> z_i;
  //  Eigen::Vector<double, D*(D+1)> xi_i;
    //double Ih = 0;
    // cout << "Into the for" << endl;
    //Eigen::VectorXd Ih(F->rows());
    Ih->setZero();

    #pragma omp parallel for
    for (int i = 0; i < F->rows(); i++) {
        Eigen::Vector<double, D*(D+1)> z_i;
        Eigen::Vector<double, D*(D+1)> xi_i;
        for (int n = 0; n < D+1; n++) {
            xi_i.segment(n*D, D) = (*Vc)((*F)(i,n), Eigen::all);
        }

        z_i = z.segment(D*(D+1)*i, D*(D+1));

        // assert(x_i.isApprox(xi_i));

        (*Ih)(i) += newtonOptSimplex(i, z_i, xi_i, 2);
        // Ih = newtonOptSimplex(i, z_i, xi_i, 10);

        z.segment(D*(D+1)*i, D*(D+1)) = z_i;

        // Fix any boundary points
        for (int n = 0; n < D+1; n++) {
            if ((*boundaryMask)((*F)(i,n))) {
	        z.segment(D*(D+1)*i+n*D, D) = x.segment((*F)(i,n)*D, D);
            }
        }
    }

    return Ih->sum();
    // cout << "Ih = " << Ih << endl;
}

template <int D>
void Mesh<D>::copyX(Eigen::VectorXd &tar) {
    int cols = Vp->cols();
    for (int i = 0; i < Vp->rows(); i++) {
        for (int j = 0; j < cols; j++) {
            tar(i*cols+j) = (*Vp)(i, j);
        }
    }
}

template <int D>
void Mesh<D>::setUp() {

    // Update the mesh in the interpolator.
    // cout << "Updating the mesh" << endl;
    mapEvaluator->updateMesh((*this->Vc), (*this->F));
    // cout << "FINISHED Updating the mesh" << endl;
    // cout << "inteprolating the monitor function" << endl;
    mapEvaluator->interpolateMonitor(*Mon);

    // cout << "FINSIHED interpolating the monitor function" << endl;

    // Mon->interpolateMonitor(*mapEvaluator, *Vc, *F, *monitorEvals);

    // // Evaluate and smooth the monitor function.
    // cout << "Eval at vertices" << endl;
    // Mon->evaluateAtVertices(*Vc, *F, *monitorEvals);;
    // cout << "FINISHED Eval at vertices" << endl;
    // cout << "smmooth monitor" << endl;
    // Mon->smoothMonitor(NUM_SMOOTH, (*this->Vc), (*this->F), *monitorEvals, *mapEvaluator);
    // cout << "Finished setting up" << endl;
}

template <int D>
void Mesh<D>::updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x) {
    int cols = Vp->cols();
    for (int i = 0; i < Vp->rows(); i++) {
        for (int j = 0; j < cols; j++) {
            (*Vp)(i, j) = x(i*cols+j);
        }
    }
}

template <int D>
void Mesh<D>::outputSimplices(const char *fname) {
    std::ofstream outFile;
    outFile.open(fname);

    for (int i = 0; i < F->rows(); i++) {
        outFile << (*F)(i, 0) << ", " << (*F)(i, 1) << ", " << (*F)(i, 2) << endl;
    }

    outFile.close();
}

template <int D>
void Mesh<D>::outputPoints(const char *fname) {
    std::ofstream outFile;
    outFile.open(fname);

    for (int i = 0; i < Vp->rows(); i++) {
        outFile << (*Vp)(i, 0) << ", " << (*Vp)(i, 1) << endl;
    }

    outFile.close();
}

template <int D>
void Mesh<D>::printDiff() {
    cout << "Difference after step: " << (*Vp - *Vc).norm() << endl;
}

template <int D>
Mesh<D>::~Mesh() {

    // delete Vc;
    // delete Vp;
    // delete F;
    delete DXpU;

    delete M;
    delete Dmat;
    delete W;

    delete I_wx;
    // delete boundaryMask;
}

// explicit instantiation for each dimension of interest
template class Mesh<2>;
template class Mesh<3>;
