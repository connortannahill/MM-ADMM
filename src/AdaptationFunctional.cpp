#include "AdaptationFunctional.h"
#include <Eigen/Dense>
#include "MonitorFunction.h"
#include <vector>
#include <iostream>
#include <Eigen/StdVector>

using namespace std;

template <int D>
AdaptationFunctional<D>::AdaptationFunctional(const AdaptationFunctional &obj) {
    Vc = obj.Vc;
    Vp = obj.Vp;
    DXpU = obj.DXpU;

    this->M = obj.M;
    // this->boundaryNode = obj.boundaryNode;
    cout << "In CTOR" << endl;
    assert(false);

    // assert(false);
    // mPre = new vector<Eigen::Matrix<double, D, D>>(D+1);
    // for (int i = 0; i < D+1; i++) {
    //     mPre->at(i) = obj.mPre->at(i);
    // }
    // for (int i = 0; i < D+1; i++) {
    //     // Eigen::Matrix<double, D, D> temp = Eigen::Matrix<double, D, D>::Zero();
    //     // mPre->push_back(temp);
    // }
    // assert(false);
}

/**
 * Bring the specified id to the front of the vector via a swap
*/
template <int D>
void swap(int id, Eigen::Vector<int, D+1> &ids) {
    for (int i = 0; i < ids.size(); i++) {
        if (ids[i] == id) {
            ids[i] = ids[0];
            ids[0] = id;
            break;
        }
    }
}

template <int D>
AdaptationFunctional<D>::AdaptationFunctional( Eigen::MatrixXd &Vc,
        Eigen::MatrixXd &Vp, Eigen::MatrixXi &F,
        Eigen::VectorXd &DXpU, MonitorFunction<D> *m, double w) {
    
    if (D == -1) {
        cout << "Must specify DIM >= 1!" << endl;
        assert(D == -1);
    }
    
    this->w = w;
    // this->boundaryNode = boundaryNode;

    this->DXpU = &DXpU;

    this->Vc = &Vc;
    this->Vp = &Vp;
    this->DXpU = &DXpU;

    this->M = m;

    // this->mPre = new vector<Eigen::Matrix<double, D, D>>(D+1);
}


/**
 * Method for computing the block gradient of a simplex for the objective function
*/
template <int D>
inline double AdaptationFunctional<D>::blockGrad(int zId, Eigen::Vector<double, D*(D+1)> &z,
            Eigen::Vector<double, D*(D+1)> &xi,
            Eigen::Vector<double, D*(D+1)> &grad,
            MeshInterpolator<D> &interp, bool computeGrad, bool regularize) {
    double Ih = 0.0;
    double detFJ;
    Eigen::Vector<double,D> gradSimplex;
    // Eigen::Vector<double,D*(D+1)> gradLocal;

    Eigen::Matrix<double,D,D> E;
    Eigen::Matrix<double,D,D> FJ;
    Eigen::Matrix<double,D,D> Ehat;
    Eigen::Matrix<double,D,D> M;
    Eigen::Vector<double,D> xK;
    Eigen::Matrix<double,D,D> vLoc;
    Eigen::Vector<double,D> xTemp;
    Eigen::Matrix<double,D,D> Einv;
    Eigen::Matrix<double,D,D> dGdJ;
    Eigen::Matrix<double,D,D> Mt0;
    Eigen::Matrix<double,D,D> Mtn;
    Eigen::Vector<double,D> dGdX;
    Eigen::Matrix<double,D,D> dGdM;
    Eigen::Vector<double,D> basisComb;
    Eigen::Vector<int, D+1> ids;

    double dFact;
    if (D == 2) {
        dFact = 2.0;
    } else if (D == 3) {
        dFact = 6.0;
    }
    double absK;

    // Build the centroid
    xK = z.segment(0, D);
    for (int n = 1; n <= D; n++) {
        xK += z.segment(D*n, D);
    }
    xK /= ((double) D + 1.0);

    // Interpolate the monitor function
    interp.evalMonitorOnGrid(xK, M);
    Eigen::Matrix<double, D, D> Minv(M.inverse());
    vector<Eigen::Matrix<double, D, D>> mPre(D+1);

    Eigen::Matrix<double, D, D> mTemp;
    for (int i = 0; i < D+1; i++) {
        xTemp = z.segment(i*D, D);
        interp.evalMonitorOnGrid(xTemp, mTemp);
        mPre.at(i) = mTemp;
    }

    double G, dGddet;

    // Compute the edge matrix
    int j = 0;
    int i = 0;
    for (int n = (i+1)%(D+1); n != i; n = (n+1)%(D+1)) {
        E.col(j)    = z.segment(D*n, D)  - z.segment(D*i, D);
        Ehat.col(j) = xi.segment(D*n, D) - xi.segment(D*i, D);
        j++;
    }

    // The element volume
    double Edet = E.determinant();

    Einv = E.inverse();

    // Approximate Jacobian
    FJ = Ehat * Einv;
    detFJ = FJ.determinant();
    double d = (double) D;
    double p = 1.5;
    double theta = 1.0/3.0;

    Eigen::Matrix<double, D, D> FJt = FJ.transpose();
    Eigen::Matrix<double, D, D> MinvJt = Minv * FJt;
    Eigen::Matrix<double, D, D> JMJt = FJ * MinvJt;
    double trJMJt = JMJt.trace();
    double detM = sqrt(1.0 / Minv.determinant());

    G = theta * detM * pow(trJMJt, d*p/2.0)
        + (1.0 - 2.0*theta) * pow(d, d*p/2.0) * detM * pow(detFJ/detM, p);
    // G = this->G(FJ, detFJ, M, xK);

    absK = abs(Edet/dFact);

    if (!computeGrad) {
        return absK * G + 0.5*w*w*( (*DXpU).segment(D*(D+1)*zId, D*(D+1)) - z ).squaredNorm();
    }

    // this->dGdJ(FJ, detFJ, M, xK, dGdJ);
    // dGddet = this->dGddet(FJ, detFJ, M, xK);
    // this->dGdM(FJ, detFJ, M, xK, dGdM);
    // this->dGdX(FJ, detFJ, M, xK, dGdX);

    dGdJ = d*p*theta*detM * pow(trJMJt, d*p/2.0 - 1) * MinvJt;
    dGddet = p*(1.0 - 2.0*theta)*pow(d, (d*p)/2.0)*pow(detM, 1.0 - p)*pow(detFJ, p-1);
    dGdM = (-0.5*theta*d*p * detM * pow(trJMJt, d*p/2.0 - 1) * Minv.transpose() * FJt * FJ * Minv)
        + (( 0.5*theta * detM * pow(trJMJt, d*p/2.0) 
        + ((0.5-theta)*(1.0-p)*pow(d, d*p/2.0)) * pow(detM, 1-p) * pow(detFJ, p) ) * Minv);
    dGdX.setZero();

    basisComb.setZero();
    j = 0;
    for (int n = (i+1)%(D+1); n != i; n = (n+1)%(D+1)) {
        basisComb += Einv.row(j) * ((dGdM * (mPre.at(n) - mPre.at(i))).trace());
        j++;
    }

    double c1 = (-G + dGddet*detFJ);
    vLoc = c1*Einv + Einv*dGdJ*FJ;
    for (int n = 0; n < D; n++) {
        vLoc.row(n) -= (basisComb)/((double) D + 1.0);
    }

    // Compute the gradient for current simplex
    gradSimplex.setZero();

    for (int n = 0; n < D; n++) {
        gradSimplex += vLoc.row(n);
    }
    gradSimplex += basisComb + dGdX;

    // Compute the gradient
    for (int l = 0; l < D; l++) {
        grad(D*i + l) = gradSimplex(l);
    }

    for (int n = 1; n < D+1; n++) {
        for (int l = 0; l < D; l++) {
            grad(D*n + l) = -vLoc(n-1, l);
        }
    }

    grad *= absK;

    // cout << "IN BLOCKGRAD" << endl;
    // cout << "z_i " << z.transpose() << endl;
    // cout << "xi_i " << z.transpose() << endl;
    // cout << "grad norm z " << grad.norm() << endl;
    // cout << "grad " << grad.transpose() << endl;
    // cout << "OUT BLOCKGRAD" << endl;

    // Update the energy
    if (i == 0) {
        Ih = absK * G;
    }

    // Now add the constraint regularization
    if (regularize) {
        Ih += 0.5*w*w*( (*DXpU).segment(D*(D+1)*zId, D*(D+1)) - z ).squaredNorm();

        grad += w*w*(-(*DXpU).segment(D*(D+1)*zId, D*(D+1)) + z);
    }

    return Ih;
}

template <int D>
AdaptationFunctional<D>::~AdaptationFunctional() {
}

template class AdaptationFunctional<2>;
template class AdaptationFunctional<3>;
