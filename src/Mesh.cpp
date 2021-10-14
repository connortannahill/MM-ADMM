#include "Mesh.h"
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
#include <unistd.h>
#include <algorithm>
#include <Eigen/LU>

#ifdef THREADS
#include <omp.h>
#endif

using namespace std;

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
void Mesh<D>::buildFaceList() {
    // Construct face connections vector
    faceConnects = new vector<set<int>>(Vp->rows());

    vector<vector<int>> faceTemp;
    // Append vertices of each boundary face to this vector
    for (int i = 0; i < F->rows(); i++) {
        // Take sum of boundary labels
        int sum = 0;
        for (int j = 0; j < D+1; j++) {
            sum += (boundaryMask->at((*F)(i, j)) != NodeType::INTERIOR) ? 1 : 0;
        }

        // If two of the virtices of this tet is on the face, it is a boundary simplex.
        if (sum == D) {
            vector<int> pnts;
            for (int j = 0; j < D+1; j++) {
                if (boundaryMask->at((*F)(i, j)) != NodeType::INTERIOR) {
                    pnts.push_back((*F)(i, j));
                }
            }

            faceTemp.push_back(pnts);
        }
    }

    // Convert this into an Eigen matrix for consistency.
    faceList = new Eigen::MatrixXi(faceTemp.size(), D);
    for (uint32_t i = 0; i < faceTemp.size(); i++) {
        for (int j = 0; j < D; j++) {
            (*faceList)(i, j) = (faceTemp.at(i)).at(j);
        }
    }

    // Build the maps for each point to keep track of boundary face membership
    for (int i = 0; i < faceList->rows(); i++) {
        for (int j = 0; j < D; j++) {
            faceConnects->at((*faceList)(i, j)).insert(i);
        }
    }
}

static int sgn(double val) {
    return (double(0) < val) - (val < double(0));
}

template <int D>
void Mesh<D>::projection2D(int nodeId, Eigen::Vector<double, D> &nodeIn) {
    if (D == 3) {
        assert(false);
    }
    Eigen::Vector2d node;
    node(0) = nodeIn(0);
    node(1) = nodeIn(1);

    Eigen::Vector2d minNodeProj = node;

    double minDist = INFINITY;
    set<int> faceIds(faceConnects->at(nodeId));

    for (auto faceId = faceIds.begin(); faceId != faceIds.end(); ++faceId) {
        Eigen::Vector<int, D> pntIds((*faceList)(*faceId, Eigen::all));
        Eigen::Vector2d x1((*Vp)(pntIds(0), Eigen::all));
        Eigen::Vector2d x2((*Vp)(pntIds(1), Eigen::all));

        // Edge vector
        Eigen::Vector2d u = x2 - x1;

        // x - x1
        Eigen::Vector2d w = node - x1;

        // Projection from first edge point to line
        Eigen::Vector2d proj = (u.dot(w) / u.dot(u)) * u;

        // Compute the distance of the point to the projection
        double dist = (proj - w).norm();

        // Compute t
        double t = proj.norm() / u.norm();

        // If the distance is the least, and the projected point is within the segment,
        // keep track.
        if (dist < minDist && sgn(u(0)) == sgn(proj(0)) && sgn(u(1)) == sgn(proj(1)) &&
                t > 0.0 && t < 1.0) {
            minDist = dist;
            minNodeProj = (1.0 - t)*x1 + t*x2;
        } else if (!(sgn(u(0)) == sgn(proj(0)) && sgn(u(1)) == sgn(proj(1)))) {
            dist = (x1 - node).norm();
            if (dist < minDist) {
                minDist = dist;
                minNodeProj = x1;
            }
        } else if (t > 1.0) {
            dist = (x2 - node).norm();
            if (dist < minDist) {
                minDist = dist;
                minNodeProj = x2;
            }
        }
    }
    for (int i = 0; i < D; i++) {
        nodeIn(i) = minNodeProj(i); 
    }
}

template <int D>
void Mesh<D>::projection3D(int nodeId, Eigen::Vector<double, D> &nodeIn) {
    if (D == 2) {
        assert(false);
    }
    Eigen::Vector3d node;
    node(0) = nodeIn(0);
    node(1) = nodeIn(1);
    node(2) = nodeIn(2);

    Eigen::Vector3d b;
    Eigen::Vector3d minNodeProj;
    bool projWorked = false;

    const double CHECK_EPS = 1e-10;

    double minDist = INFINITY;
    set<int> faceIds(faceConnects->at(nodeId));

    for (auto nodePnt = faceIds.begin(); nodePnt != faceIds.end(); ++nodePnt) {
        Eigen::Vector<int, D> pntIds((*faceList)(*nodePnt, Eigen::all));
        Eigen::Vector3d pnt0((*Vp)(pntIds(0), Eigen::all));
        Eigen::Vector3d pnt1((*Vp)(pntIds(1), Eigen::all));
        Eigen::Vector3d pnt2((*Vp)(pntIds(2), Eigen::all));

        Eigen::Vector3d q = pnt0;
        Eigen::Vector3d u = pnt1 - q;
        Eigen::Vector3d v = pnt2 - q;
        Eigen::Vector3d n = (u.cross(v));

        double temp = 1.0 / n.dot(n);
        Eigen::Vector3d w = node - q;

        b(2) = (u.cross(w)).dot(n) * temp;
        b(1) = (w.cross(v)).dot(n) * temp;
        b(0) = 1.0 - b(1) - b(2);

        Eigen::Vector3d nodeProj = b(0)*(*Vp)(pntIds(0), Eigen::all)
            + b(1)*(*Vp)(pntIds(1), Eigen::all)
            + b(2)*(*Vp)(pntIds(2), Eigen::all);

        double dist = (nodeProj - node).norm();
        if (dist < minDist && b(0) >= CHECK_EPS && b(1) >= CHECK_EPS && b(2) >= CHECK_EPS) {
            minDist = dist;
            minNodeProj = nodeProj;
            projWorked = true;
        }
    }

    // If the projection is on the edge, good. Otherwise, the node does not move.
    if (projWorked) {
        for (int i = 0; i < D; i++) {
            nodeIn(i) = minNodeProj(i); 
        }
    }
}

template <int D>
void Mesh<D>::projectOntoBoundary(int nodeId, Eigen::Vector<double, D> &node) {
    if (D == 2) {
        projection2D(nodeId, node);
    } else if (D == 3) {
        projection3D(nodeId, node);
    }
}

template <int D>
void Mesh<D>::reOrientElements(Eigen::MatrixXd &Xp, Eigen::MatrixXi &F) {
    Eigen::Matrix<double,D,D> E;

    // Compute the edge matrix
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < D; j++) {
            E.col(j)    = Xp(F(i, j+1), Eigen::all)  - Xp(F(i, 0), Eigen::all);
        }

        if (E.determinant() < 0) {
            int temp = F(i, 2);
            F(i, 2) = F(i, 1);
            F(i, 1) = temp;
        }
    }

}

template <int D>
void Mesh<D>::meshInit(Eigen::MatrixXd &Xc, Eigen::MatrixXd &Xp, 
            Eigen::MatrixXi &F, vector<NodeType> &boundaryMask,
            MonitorFunction<D> *Mon, int numThreads, double rho, double tau) {

    if (compMesh) {
        this->Vc = &Xc;
    } else {
        this->Vc = nullptr;
    }

    this->Vp = &Xp;

    // Should only use positive orientation elements
    reOrientElements(Xp, F);

    this->F = &F;
    this->boundaryMask = &boundaryMask;

    this->tau = tau;

    this->Ih = new Eigen::VectorXd(F.rows());

    buildFaceList();

#ifdef THREADS
    omp_set_num_threads(numThreads);
#endif

    // Create mesh interpolator
    mapEvaluator = new MeshInterpolator<D>();
    mapEvaluator->updateMesh((*this->Vp), (*this->F));
    mapEvaluator->interpolateMonitor(*Mon);
    this->nPnts = Xp.rows();

    this->Mon = Mon;

    // Trivial edge matrix
    int nPnts = Vp->rows();

    // Build the monitor simplex values
    monitorEvals = new Eigen::MatrixXd(Xp.rows(), D*D);

    // Build the mass matrix
    Eigen::VectorXd m(Eigen::VectorXd::Constant(nPnts*D, tau));
    this->m = tau;
    buildMassMatrix(m);

    buildDMatrix();
    this->w = sqrt(rho);
    buildWMatrix(this->w);

    DXpU = new Eigen::VectorXd((*Dmat)*(m));

    hessInvs = new vector<Eigen::Matrix<double, D*(D+1), D*(D+1)>>(this->F->rows());
    gradCurrs = new vector<Eigen::Vector<double, D*(D+1)>>(this->F->rows());

    // Set the Hessians to the identity intiially
    Eigen::Matrix<double, D*(D+1), D*(D+1)> eye
        = Eigen::Matrix<double, D*(D+1), D*(D+1)>::Identity(D*(D+1), D*(D+1));
    for (uint32_t i = 0; i < hessInvs->size(); i++) {
        hessInvs->at(i) = eye;
    }

    // Create the functional for each vertex
    if (compMesh) {
        I_wx = new HuangFunctional<D>(*Vc, *Vp, *(this->F), *DXpU, Mon, this->w, 0.0, 0.0);
    } else {
        I_wx = new HuangFunctional<D>(*Vp, *(this->F), *DXpU, Mon, this->w, 0.0, 0.0);
    }
}

template <int D>
Mesh<D>::Mesh(Eigen::MatrixXd &Xp, Eigen::MatrixXi &F, vector<NodeType> &boundaryMask,
            MonitorFunction<D> *Mon, int numThreads, double rho, double tau) {
    
    this->compMesh = false;

    Eigen::MatrixXd Xc;
    meshInit(Xc, Xp, F, boundaryMask, Mon, numThreads, rho, tau);

}

template <int D>
Mesh<D>::Mesh(Eigen::MatrixXd &Xc, Eigen::MatrixXd &Xp, Eigen::MatrixXi &F, vector<NodeType> &boundaryMask,
            MonitorFunction<D> *Mon, int numThreads, double rho, double tau) {
    
    compMesh = true;

    meshInit(Xc, Xp, F, boundaryMask, Mon, numThreads, rho, tau);
}

template <int D>
double Mesh<D>::computeEnergy(Eigen::VectorXd &x) {
    double Ih = 0.0;

    Eigen::Vector<double, D*(D+1)> x_i;
    Eigen::Vector<double, D*(D+1)> xi_i;
    Eigen::Vector<double, D*(D+1)> gradSimp;
    Eigen::Vector<double, D> pnt;

#ifdef THREADS
    #pragma omp parallel for
#endif
    for (int i = 0; i < F->rows(); i++) {
        if (compMesh) {
            for (int n = 0; n < D+1; n++) {
                pnt = (*Vc)((*F)(i,n), Eigen::all);
                
                for (int l = 0; l < D; l++) {
                    xi_i(n*D+l) = pnt(l);
                }
            }
        }

        for (int n = 0; n < D+1; n++) {
            pnt = (*Vp)((*F)(i,n), Eigen::all);

            for (int l = 0; l < D; l++) {
                x_i(n*D+l) = pnt(l);
            }
        }

        // Compute the local gradient
        Ih += I_wx->blockGrad(i, x_i, xi_i, gradSimp, *mapEvaluator, false, false);
    }

    return Ih;
}

template <int D>
double Mesh<D>::eulerStep(Eigen::VectorXd &x, Eigen::VectorXd &grad) {
    Eigen::Vector<double, D*(D+1)> x_i;
    Eigen::Vector<double, D*(D+1)> xi_i;
    Eigen::Vector<double, D*(D+1)> gradSimp;
    Eigen::Vector<double, D> pnt;

    grad.setZero();
    double Ihorig = 0.0;

    for (int i = 0; i < F->rows(); i++) {
        if (compMesh) {
            for (int n = 0; n < D+1; n++) {
                pnt = (*Vc)((*F)(i,n), Eigen::all);
                
                for (int l = 0; l < D; l++) {
                    xi_i(n*D+l) = pnt(l);
                }
            }
        }

        for (int n = 0; n < D+1; n++) {
            pnt = x.segment(D*(*F)(i,n), D);

            for (int l = 0; l < D; l++) {
                x_i(n*D+l) = pnt(l);
            }
        }

        // Compute the local gradient
        Ihorig += I_wx->blockGrad(i, x_i, xi_i, gradSimp, *mapEvaluator, true, false);

        // For place into the respective gradients
        for (int n = 0; n < D+1; n++) {
            int off = (*F)(i,n);
                
            if (boundaryMask->at((*F)(i,n)) == NodeType::INTERIOR) {
                grad.segment(D*off, D) += gradSimp.segment(D*n, D);
            }
        }
    }

    cout << "in Euler Ihorig = " << Ihorig << endl;
    // assert(false);

    return Ihorig;
}

// TODO: add external forces (would be handled in this object, as the meshing object should control its own props)
template <int D>
double Mesh<D>::predictX(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar) {
    double dtNew = 0;
    if (!stepTaken) {
        // If the step is not taken we do Euler's method
        Eigen::VectorXd grad(D*Vp->rows());

        // Compute the initial guess
        double Ihorig = eulerStep(x, grad);

        // assert(false);
        // cout << "dt = " << dt << endl;
        // cout << "tau = " << tau << endl;
        // cout << "grad = " << grad.transpose() << endl;
        // assert(false);
        xBar = x - (dt / tau) * grad;

        // Compute energy at the current time level along with the updated gradient
        double Ih = eulerStep(xBar, grad);
        // assert(false);
        double rate = (Ih - Ihorig)/Ihorig;

        dtNew = dt;
        while (abs(rate) > 1e-2 && rate > 0) {
            dtNew /= 2.0;
            xBar = x - (dtNew / tau) * grad;
            Ih = eulerStep(xBar, grad);
            rate = (Ih - Ihorig)/Ihorig;
        }
    } else {
        // Compute the gradient at each point
        xBar = 2.0*x - xPrev;
        dtNew = dt;
    }

    return dtNew;
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
    // Temp matrix (we insert the transpose) and we keep it this way for storage reasons
    Eigen::SparseMatrix<double> Dt(D*Vp->rows(), D*(D+1)*F->rows());


    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(D*(D+1)*F->rows());

    // Add the spring constraints to the system
    int pntIds[D+1];
    int colStarts[D+1];
    int rowStarts[D+1];
    // Note: the is columns in the matrix transposed
    for (int col = 0; col < F->rows(); col++) {
        // Extract the ids
        for (int i = 0; i < D+1; i++)  {
            pntIds[i] = (*F)(col, i);
        }

        colStarts[0] = D*(D+1)*col;
        for (int i = 1; i < D+1; i++) {
            colStarts[i] = colStarts[i-1]+D;
        }

        for (int i = 0; i < D+1; i++) {
            rowStarts[i] = pntIds[i]*D;
        }

        // First block
        for (int n = 0; n < D; n++) {
            for (int m = 0; m < D; m++) {
                if (n == m) {
                    for (int i = 0; i < D+1; i++) {
                        tripletList.push_back(T(n+rowStarts[i], m+colStarts[i], 1.0));
                    }
                } else {
                    for (int i = 0; i < D+1; i++) {
                        tripletList.push_back(T(n+rowStarts[i], m+colStarts[i], 0.0));
                    }
                }
            }
        }
    }

    Dt.setFromTriplets(tripletList.begin(), tripletList.end());

    // Create the D matrix
    Dmat = new Eigen::SparseMatrix<double>(Dt.transpose());
}

/**
 * BFGS quasi-Newton's method over a single simplex
*/
template <int D>
inline double Mesh<D>::bfgsOptSimplex(int zId, Eigen::Vector<double, D*(D+1)> &z,
        Eigen::Vector<double, D*(D+1)> &xi, int nIter, double tol, bool firstStep) {
    double h = 2.0*sqrt(std::numeric_limits<double>::epsilon());

    Eigen::Vector<double, D*(D+1)> zPurt;
    Eigen::Vector<double, D*(D+1)> Gkp1;
    Eigen::Vector<double, D*(D+1)> Gk;

    double Ix = 0;

    // If this is the first step for this simplex, compute the initial
    // Hessian and gradient
    Ix = I_wx->blockGrad(zId, z, xi, Gk, *mapEvaluator, true, true);
    if (!hessComputed) {
        zPurt = z;
        for (int i = 0; i < D*(D+1); i++) {
            // Compute purturbation
            zPurt(i) += h;

            // Compute gradient at purturbed point
            I_wx->blockGrad(zId, zPurt, xi, Gkp1, *mapEvaluator, true, true);
            hessInvs->at(zId).col(i) = (Gkp1 - Gk)/h;

            zPurt(i) = z(i);
        }
        hessInvs->at(zId) = hessInvs->at(zId).inverse();
    }

    if (Gk.norm() < tol) {
        return Ix;
    }

    // Get the existing information
    Eigen::Matrix<double, D*(D+1), D*(D+1)> Bkinv = hessInvs->at(zId);
    Eigen::Vector<double, D*(D+1)> pk;
    Eigen::Vector<double, D*(D+1)> yk;
    Eigen::Vector<double, D*(D+1)> zTemp;
    double c1, c2;

    int iter;
    for (iter = 0; iter < nIter; iter++) {
        // Compute the Newton direction
        pk = - Bkinv * Gk; 

        // Update the solution
        z += pk;

        // Update the secant direction, first computing the new gradient
        Ix = I_wx->blockGrad(zId, z, xi, Gkp1, *mapEvaluator, true, true);

        if (Gkp1.norm() < tol) {
            Gk = Gkp1;
            break;
        }

        yk = Gkp1 - Gk; 

        // Update information for the next step
        c2 = pk.dot(yk);
        c1 = (pk.dot(yk) + yk.dot(Bkinv*yk))/(pow(c2, 2.0));
        Bkinv += c1 * (pk * pk.transpose()) - Bkinv*(yk*pk.transpose()) / c2 - pk*(yk.transpose()*Bkinv) / c2; 
        Gk = Gkp1;
    }

    hessInvs->at(zId) = Bkinv;

    // Return the updated objective function value on this simplex
    return Ix;
}

/**
 * Newton's method over a single simplex
*/
template <int D>
double Mesh<D>::newtonOptSimplex(int zId, Eigen::Vector<double, D*(D+1)> &z,
        Eigen::Vector<double, D*(D+1)> &xi, int nIter, double tol) {
    double h = 2.0*sqrt(std::numeric_limits<double>::epsilon());

    Eigen::Vector<double, D*(D+1)> zPurt;
    Eigen::Vector<double, D*(D+1)> gradZ;
    Eigen::Vector<double, D*(D+1)> gradZPurt;
    Eigen::Matrix<double, D*(D+1), D*(D+1)> hess;
    Eigen::Matrix<double, D*(D+1), D*(D+1)> hessInv;
    Eigen::Vector<double, D*(D+1)> p;
    Eigen::Vector<double, D*(D+1)> gradTemp;

    double Ix = 0;
    double IxPrev = INFINITY;

    for (int iters = 0; iters < nIter; iters++) {
        Ix = I_wx->blockGrad(zId, z, xi, gradZ, *mapEvaluator, true, true);

        if (iters != 0 && (abs((IxPrev - Ix) / IxPrev) < tol)) {
            break;
        }

        zPurt = z;

        // Compute the Hessian column-wise (we use fixed Hessians for all
        // iterations other than the first)
        if (iters == 0) {
            for (int i = 0; i < D*(D+1); i++) {
                // Compute purturbation
                zPurt(i) += h;

                // Compute gradient at purturbed point
                I_wx->blockGrad(zId, zPurt, xi, gradZPurt, *mapEvaluator, true, true);
                hess.col(i) = (gradZPurt - gradZ)/h;

                zPurt(i) = z(i);
            }
        }


        // Compute the Newton direction
        p = hess.lu().solve(-gradZ);
        z += p;

        IxPrev = Ix;
    }

    return Ix;

}

template <int D>
double Mesh<D>::prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DXpU, Eigen::VectorXd &z) {
    // Copy DXpU address to local pointer
    *this->DXpU = DXpU;

    // Run Newton's method on each simplex
    Ih->setZero();

#ifndef THREADS
    Eigen::Vector<double, D*(D+1)> z_i;
    Eigen::Vector<double, D*(D+1)> xi_i;
    Eigen::Vector<double, D> pnt;
    Eigen::Vector<double, D> zTemp;
#endif

#ifdef THREADS
    #pragma omp parallel for
#endif
    for (int i = 0; i < F->rows(); i++) {

#ifdef THREADS
        Eigen::Vector<double, D*(D+1)> z_i;
        Eigen::Vector<double, D*(D+1)> xi_i;
        Eigen::Vector<double, D> pnt;
        Eigen::Vector<double, D> zTemp;
#endif
        if (compMesh) {
            for (int n = 0; n < D+1; n++) {
                pnt = (*Vc)((*F)(i,n), Eigen::all);
                
                for (int l = 0; l < D; l++) {
                    xi_i(n*D+l) = pnt(l);
                }
            }
        }

        z_i = z.segment(D*(D+1)*i, D*(D+1));

        (*Ih)(i) = bfgsOptSimplex(i, z_i, xi_i, 100, 1e-6, hessComputed);

        for (int l = 0; l < D*(D+1); l++) {
            z(D*(D+1)*i+l) = z_i(l);
        }

        // Moving boundary points
        for (int n = 0; n < D+1; n++) {
            if (boundaryMask->at((*F)(i,n)) == NodeType::BOUNDARY_FREE) {
                zTemp = z.segment(D*(D+1)*i+n*D, D);
                projectOntoBoundary((*F)(i,n), zTemp);
                for (int l = 0; l < D; l++) {
                    z(D*(D+1)*i+n*D + l) = zTemp(l);
                }
            } else if (boundaryMask->at((*F)(i,n)) == NodeType::BOUNDARY_FIXED) { 
                z.segment(D*(D+1)*i+n*D, D) = x.segment((*F)(i,n)*D, D);
            }
        }

    }
    hessComputed = true;

    return Ih->sum();
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
    if (!stepTaken) {
        // Update the mesh in the interpolator.
        mapEvaluator->updateMesh((*this->Vp), (*this->F));
        mapEvaluator->interpolateMonitor(*Mon);
        // stepTaken = true;
    }
}

template <int D>
void Mesh<D>::updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x) {
    int cols = Vp->cols();
    for (int i = 0; i < Vp->rows(); i++) {
        for (int j = 0; j < cols; j++) {
            // (*Vc)(i, j) = xPrev(i*cols+j);
            (*Vp)(i, j) = x(i*cols+j);
        }
    }
}

template <int D>
void Mesh<D>::outputBoundaryNodes(const char *fname) {
    std::ofstream outFile;
    outFile.open(fname);

    for (int i = 0; i < Vp->rows(); i++) {
        if (boundaryMask->at(i) != INTERIOR) {
            for (int j = 0; j < D-1; j++) {
                outFile << (*Vp)(i, j) << ", ";
            }
            outFile << (*Vp)(i, D-1) << endl;

        }
    }

    outFile.close();
}

template <int D>
void Mesh<D>::outputSimplices(const char *fname) {
    std::ofstream outFile;
    outFile.open(fname);

    for (int i = 0; i < F->rows(); i++) {
        for (int j = 0; j < D; j++) {
            outFile << (*F)(i, j) << ", ";
        }
        outFile << (*F)(i, D) << endl;
    }

    outFile.close();
}

template <int D>
void Mesh<D>::outputPoints(const char *fname) {
    std::ofstream outFile;
    outFile.open(fname);

    for (int i = 0; i < Vp->rows(); i++) {
        for (int j = 0; j < D-1; j++) {
            outFile << (*Vp)(i, j) << ", ";
        }
        outFile << (*Vp)(i, D-1) << endl;
    }

    outFile.close();
}

template <int D>
void Mesh<D>::printDiff() {
    cout << "Difference after step: " << (*Vp - *Vc).norm() << endl;
}

template <int D>
Mesh<D>::~Mesh() {

    delete DXpU;
    delete Ih;

    delete M;
    delete Dmat;
    delete W;

    delete faceList;
    delete faceConnects;

    delete I_wx;
    delete hessInvs;
    delete gradCurrs;
}

// explicit instantiation for each dimension of interest
template class Mesh<2>;
template class Mesh<3>;
