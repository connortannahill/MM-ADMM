#include "Mesh2D.h"
#include "Assembly.h"
#include <Eigen/Dense>
#include <unordered_map>
#include "MonitorFunction.h"
#include <iostream>
#include <LBFGS.h>
#include <fstream>
#define DIM 2
// #include "igl/triangle/triangulate.h"

using namespace std;
using namespace LBFGSpp;

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

// TODO: another constructor for arb. domains should take in a traingulation from the user.
Mesh2D::Mesh2D(MonitorFunction *M, unordered_map<string, double> params) : Assembly(params) {
    this->setD(2);

    this->nx = (int) params["nx"];
    this->ny = (int) params["ny"];

    this->xa = params["xa"];
    this->xb = params["xb"];
    this->ya = params["ya"];
    this->yb = params["yb"];

    this->rho = params["rho"];
    this->tau = params["tau"];

    this->hx = (xb - xa)/((double)this->nx);
    this->hy = (yb - ya)/((double)this->ny);

    this->M = M;

    // Build up the triangulations

    // Trivial edge matrix
    int nPnts = (nx+1)*(ny+1) + nx*ny;
    this->setNPnts(nPnts);
    Vc = new Eigen::MatrixXd(this->nPnts, d);
    Vp = new Eigen::MatrixXd(this->nPnts, d);
    F = new Eigen::MatrixXi(4*nx*ny, d+1);

    int off = 0;
    for (int j = 0; j <= ny; j++) {
        for (int i = 0; i <= nx; i++) {
            (*Vc)(off, 0) = hx*i;
            (*Vc)(off, 1) = hy*j;

            (*Vp)(off, 0) = hx*i;
            (*Vp)(off, 1) = hy*j;

            off++;
        }
    }

    // Append the midoints
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            (*Vc)(off, 0) = hx*i + hx/2.0;
            (*Vc)(off, 1) = hy*j + hy/2.0;

            (*Vp)(off, 0) = hx*i + hx/2.0;
            (*Vp)(off, 1) = hy*j + hy/2.0;

            off++;

        }
    }

    int stride = (nx+1) * (ny+1);

    off = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            // Left
            (*F)(off, 0) = i            + j*(nx+1);
            (*F)(off, 1) = i            + (j+1)*(nx+1);
            (*F)(off, 2) = stride + i   + j*nx;
            // (*F)(off, 2) = i+1 + (j+1)*(nx+1);

            off++;

            // Top
            (*F)(off, 0) = i            + (j+1)*(nx+1);
            (*F)(off, 1) = i+1          + (j+1)*(nx+1);
            (*F)(off, 2) = stride + i   + j*nx;
            // (*F)(off, 2) = i+1 + (j+1)*(nx+1);

            off++;

            // Right
            (*F)(off, 0) = i+1          + (j+1)*(nx+1);
            (*F)(off, 1) = i+1          + j*(nx+1);
            (*F)(off, 2) = stride + i   + j*nx;
            // (*F)(off, 2) = i+1 + (j+1)*(nx+1);

            off++;

            // Bot
            (*F)(off, 0) = i            + j*(nx+1);
            (*F)(off, 1) = i+1          + j*(nx+1);
            (*F)(off, 2) = stride + i   + j*nx;
            // (*F)(off, 2) = i+1 + (j+1)*(nx+1);

            off++;

            // (*F)(off, 0) = i   + j*(nx+1);
            // (*F)(off, 1) = i+1 + j*(nx+1);
            // (*F)(off, 2) = i+1 + (j+1)*(nx+1);
            // off++;
        }
    }

    // Build the mass matrix
    Eigen::VectorXd m(Eigen::VectorXd::Constant(nPnts*d, tau));
    buildMassMatrix(m);

    buildDMatrix();
    this->w = sqrt(rho);
    // double w = sqrt(rho);
    buildWMatrix(w);

    DXpU = new Eigen::VectorXd((*D)*(m));

    // Create the functional for each vertex
    // I_wx = new vector<HuangFunctional>();
    I_wx = new HuangFunctional(DIM, false, *Vc, *Vp, *F, *DXpU, M, w, 0.0, 0.0);

    boundaryMask = new Eigen::VectorXi(d*Vp->rows());
    boundaryMask->setZero();

    for (int i = 0; i < (nx+1)*(ny+1); i++) {
        int iOff = i % (nx+1);
        int jOff = i / (ny+1);
        bool boundaryPnt = (iOff == 0) || (iOff == nx) || (jOff == 0) || (jOff == ny);

        if (boundaryPnt) {
            (*boundaryMask)(i) = 1;
        }

        // I_wx->push_back(HuangFunctional(d, i, boundaryPnt, *Vc, *Vp, *F, *DXpU, M,
        //                                 w, params["theta"],  params["p"]));
    }

    // // Add the stride points
    // for (int i = (nx+1)*(ny+1); i < nPnts; i++) {
    //     I_wx->push_back(HuangFunctional(d, i, false, *Vc, *Vp, *F, *DXpU, M,
    //                                     w, params["theta"],  params["p"]));
    // }

    // assert(I_wx->size() == nPnts);

    // cout << "Finished setting up!" << endl;
}

// TODO: add external forces (would be handled in this object, as the meshing object should control its own props)
void Mesh2D::predictX(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar) {
    xBar = 2.0*x - xPrev;
}

/**
 * Default implementation, W is the identity matrix
*/
void Mesh2D::buildWMatrix(double w) {
    this->w = w;

    if (!wAlloc) {
        W = new Eigen::SparseMatrix<double>(F->rows()*d*(d+1), F->rows()*d*(d+1));
        wAlloc = true;
    } else {
        assert(false);
    }

    W->reserve(Eigen::VectorXd::Constant(F->rows()*d*(d+1), 1));

    for (int i = 0; i < W->rows(); i++) {
        W->insert(i, i) = w;
    }
}

void Mesh2D::buildDMatrix() {
    if (dAlloc) {
        cout << "D already set!!" << endl;
        assert(false);
    }

    // Temp matrix (we insert the transpose) and we keep it this way for storage reasons
    Eigen::SparseMatrix<double> Dt(d*Vp->rows(), d*(d+1)*F->rows());

    // Reserve the memory
    Dt.reserve(Eigen::VectorXi::Constant(d*(d+1)*F->rows(), 1));

    // Add the spring constraints to the system
    int pntIds[3];
    int rowStart0, rowStart1, rowStart2, colStart0, colStart1, colStart2;
    // Note: the is columns in the matrix transposed
    for (int col = 0; col < F->rows(); col++) {
        // Extract the ids
        pntIds[0] = (*F)(col, 0);
        pntIds[1] = (*F)(col, 1);
        pntIds[2] = (*F)(col, 2);

        colStart0 = d*(d+1)*col;
        colStart1 = colStart0+d;
        colStart2 = colStart1+d;

        rowStart0 = pntIds[0]*d;
        rowStart1 = pntIds[1]*d;
        rowStart2 = pntIds[2]*d;

        // First block
        for (int n = 0; n < d; n++) {
            for (int m = 0; m < d; m++) {
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
    D = new Eigen::SparseMatrix<double>(Dt.transpose());

    dAlloc = true;
}

/**
 * Newton's method over a single simplex
*/
double Mesh2D::newtonOptSimplex(int zId, Eigen::Vector<double, DIM*(DIM+1)> &z,
        Eigen::Vector<double, DIM*(DIM+1)> &xi, int nIter) {
    double h = 2.0*sqrt(std::numeric_limits<double>::epsilon());
    const int MAX_LS = 10;

    Eigen::Vector<double, DIM*(DIM+1)> zPurt;
    Eigen::Vector<double, DIM*(DIM+1)> gradZ;
    Eigen::Vector<double, DIM*(DIM+1)> gradZPurt;
    Eigen::Matrix<double, DIM*(DIM+1), DIM*(DIM+1)> hess;
    Eigen::Vector<double, DIM*(DIM+1)> p;
    Eigen::Vector<double, DIM*(DIM+1)> gradTemp;

    double Ix;
    double Ipurt;

    for (int iters = 0; iters < nIter; iters++) {
        Ix = I_wx->blockGrad(zId, z, xi, gradZ);

        zPurt = z;

        // Compute the Hessian column-wise
        for (int i = 0; i < d*(d+1); i++) {
            // Compute purturbation
            zPurt(i) += h;

            // Compute gradient at purturbed point
            Ipurt = I_wx->blockGrad(zId, zPurt, xi, gradZPurt);

            hess.col(i) = (gradZPurt - gradZ)/h;

            zPurt(i) = z(i);
        }

        // Compute the Newton direction
        p = hess.colPivHouseholderQr().solve(-gradZ);
        zPurt = z + p;

        // Perform backtracking line search in the Hessian direction (should work if Hess is pos def)
        int lsIters = 0;
        double alpha = 1.0;
        double c1 = 0.0;
        double c2 = 0.9;
        Ipurt = I_wx->blockGrad(zId, zPurt, xi, gradTemp);

        while (Ipurt >= Ix - c1*alpha*(gradZ.dot(p)) && -p.dot(gradTemp) >= -c2*p.dot(gradZ)  && lsIters < MAX_LS) {
            alpha /= 10.0;

            zPurt = z + alpha*p;
            Ipurt = I_wx->blockGrad(zId, zPurt, xi, gradTemp);

            lsIters++;
        }

        if (lsIters == MAX_LS) {
            break;
        } else {
            z = zPurt;
        }
    }

    return Ipurt;

}

void Mesh2D::prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DXpU, Eigen::VectorXd &z) {
    // Copy DXpU address to local pointer
    *this->DXpU = DXpU;

    // Run Newton's method on each simplex
    Eigen::Vector<double, DIM*(DIM+1)> z_i;
    Eigen::Vector<double, DIM*(DIM+1)> xi_i;
    double Ih;
    for (int i = 0; i < F->rows(); i++) {
        for (int n = 0; n < d+1; n++) {
            xi_i.segment(n*d, d) = (*Vc)((*F)(i,n), Eigen::all);
        }

        z_i = z.segment(d*(d+1)*i, d*(d+1));

        // assert(x_i.isApprox(xi_i));

        Ih = newtonOptSimplex(i, z_i, xi_i, 2);
        // Ih = newtonOptSimplex(i, z_i, xi_i, 10);

        z.segment(d*(d+1)*i, d*(d+1)) = z_i;

        // assert(boundaryMask->sum() == 40);
        // cout << "Sum of boundary mask = " << boundaryMask->sum() << endl;


        // Fix any boundary points
        for (int n = 0; n < d+1; n++) {
            if ((*boundaryMask)((*F)(i,n))) {
                z.segment(d*(d+1)*i+n*d, d) = x.segment((*F)(i,n)*d, d);
            }
        }
    }
}

void Mesh2D::copyX(Eigen::VectorXd &tar) {
    int cols = Vp->cols();
    for (int i = 0; i < Vp->rows(); i++) {
        for (int j = 0; j < cols; j++) {
            tar(i*cols+j) = (*Vp)(i, j);
        }
    }
}

void Mesh2D::updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x) {
    int cols = Vp->cols();
    for (int i = 0; i < Vp->rows(); i++) {
        for (int j = 0; j < cols; j++) {
            (*Vp)(i, j) = x(i*cols+j);
        }
    }
    cout << "Difference after step: " << (*Vp - *Vc).norm() << endl;
}

void Mesh2D::outputTriangles(const char *fname) {
    std::ofstream outFile;
    outFile.open(fname);

    for (int i = 0; i < F->rows(); i++) {
        outFile << (*F)(i, 0) << ", " << (*F)(i, 1) << ", " << (*F)(i, 2) << endl;
    }

    outFile.close();
}

void Mesh2D::outputPoints(const char *fname) {
    std::ofstream outFile;
    outFile.open(fname);

    for (int i = 0; i < Vp->rows(); i++) {
        outFile << (*Vp)(i, 0) << ", " << (*Vp)(i, 1) << endl;
    }

    outFile.close();
}

void Mesh2D::printDiff() {
    cout << "Difference after step: " << (*Vp - *Vc).norm() << endl;
}

Mesh2D::~Mesh2D() {

    delete Vc;
    delete Vp;
    delete F;
    delete DXpU;

    delete I_wx;
    delete boundaryMask;
}
