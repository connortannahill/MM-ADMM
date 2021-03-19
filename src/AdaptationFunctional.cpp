#include "AdaptationFunctional.h"
#include <Eigen/Dense>
#include "MonitorFunction.h"
#include <vector>
#include <iostream>

using namespace std;

template <int D>
AdaptationFunctional<D>::AdaptationFunctional(const AdaptationFunctional &obj) {
    Vc = obj.Vc;
    Vp = obj.Vp;
    DXpU = obj.DXpU;

    this->M = obj.M;
    this->boundaryNode = obj.boundaryNode;
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
    this->boundaryNode = boundaryNode;

    this->DXpU = &DXpU;

    this->Vc = &Vc;
    this->Vp = &Vp;
    this->DXpU = &DXpU;

    this->M = m;
}

/**
 * Method for computing the block gradient of a simplex for the objective function
*/
template <int D>
double AdaptationFunctional<D>::blockGrad(int zId, Eigen::Vector<double, D*(D+1)> &z,
            Eigen::Vector<double, D*(D+1)> &xi,
            Eigen::Vector<double, D*(D+1)> &grad,
            MeshInterpolator<D> &interp) {
    double Ih = 0.0;;
    double detFJ;
    Eigen::Vector<double,D> gradSimplex;
    gradSimplex.setZero();
    grad.setZero();

    Eigen::Matrix<double,D,D> E(Eigen::Matrix<double,D,D>::Constant(0.0));
    Eigen::Matrix<double,D,D> FJ(Eigen::Matrix<double,D,D>::Constant(0.0));
    Eigen::Matrix<double,D,D> Ehat(Eigen::Matrix<double,D,D>::Constant(0.0));
    Eigen::Matrix<double,D,D> M(Eigen::Matrix<double,D,D>::Constant(0.0));
    Eigen::Vector<double,D> xK(Eigen::Vector<double,D>::Constant(0.0));
    Eigen::Matrix<double,D,D> vLoc(Eigen::Matrix<double,D,D>::Constant(0.0));
    // Eigen::Matrix<double,D,D> Mtemp(Eigen::Matrix<double,D,D>::Constant(0.0));
    Eigen::Vector<double,D> xTemp(Eigen::Vector<double,D>::Constant(0.0));
    Eigen::Matrix<double,D,D> Einv(Eigen::Matrix<double,D,D>::Constant(0.0));
    Eigen::Matrix<double,D,D> dGdJ(Eigen::Matrix<double,D,D>::Constant(0.0));
    Eigen::Matrix<double,D,D> Mt0(Eigen::Matrix<double,D,D>::Constant(0.0));
    Eigen::Matrix<double,D,D> Mtn(Eigen::Matrix<double,D,D>::Constant(0.0));
    Eigen::Vector<double,D> dGdX(Eigen::Vector<double,D>::Constant(0.0));
    Eigen::Matrix<double,D,D> dGdM(Eigen::Matrix<double,D,D>::Constant(0.0));
    Eigen::Vector<double,D> basisComb(Eigen::Vector<double,D>::Constant(0.0));
    Eigen::Vector<int, D+1> ids;

    double dFact = (double) 2.0;
    double absK;
    int d = D;

    // Build the centroid
    xK = z.segment(0, d);
    for (int n = 1; n <= d; n++) {
        xK += z.segment(d*n, d);
    }
    xK /= ((double) d + 1.0);

    // Interpolate the monitor function
    // (*this->M)(xK, M);
    interp.evalMonitorOnSimplex(zId, xK, M);
    // interp.evalMonitorAtPoint(xK, M);

    for (int i = 0; i < d+1; i++) {
        // Compute the edge matrix
        int j = 0;
        for (int n = (i+1)%(d+1); n != i; n = (n+1)%(d+1)) {
            E.col(j)    = z.segment(d*n, d)  - z.segment(d*i, d);
            Ehat.col(j) = xi.segment(d*n, d) - xi.segment(d*i, d);
            j++;
        }

        // The element volume
        double Edet = E.determinant();
        double Ehatdet = Ehat.determinant();

        absK = abs(Edet/dFact);

        Einv = E.inverse();

        // Approximate Jacobian
        FJ = Ehat * Einv;
        detFJ = Ehatdet/Edet;
        
        // Compute the gradient and Ih
        double G = this->G(FJ, detFJ, M, xK);
        this->dGdJ(FJ, detFJ, M, xK, dGdJ);
        double dGddet = this->dGddet(FJ, detFJ, M, xK);
        this->dGdX(FJ, detFJ, M, xK, dGdX);
        this->dGdM(FJ, detFJ, M, xK, dGdM);

        // Form the linear combination of the linear basis derivatives
        xTemp = z.segment(i*d, d);

        // (*this->M)(xTemp, Mt0);
        interp.evalMonitorOnSimplex(zId, xTemp, Mt0);
        // interp.evalMonitorAtPoint(xTemp, Mt0);
        
        basisComb.setZero();
        j = 0;
        for (int n = (i+1)%(d+1); n != i; n = (n+1)%(d+1)) {
            xTemp = z.segment(n*d, d);
            // interp.evalMonitorAtPoint(xTemp, Mtn);
            interp.evalMonitorOnSimplex(zId, xTemp, Mtn);
            // (*this->M)(xTemp, Mtn);
            basisComb += Einv.row(j) * ((dGdM * (Mtn - Mt0)).trace());
            j++;
        }

        vLoc = -G*Einv + Einv*dGdJ*Ehat*Einv + dGddet*detFJ*Einv;
        for (int n = 0; n < d; n++) {
            vLoc.row(n) -= (basisComb + dGdX)/((double) d + 1.0);
        }

        // Compute the gradient for current simplex
        gradSimplex.setZero();
        for (int n = 0; n < d; n++) {
            gradSimplex += vLoc.row(n);
        }
        gradSimplex += basisComb + dGdX;

        gradSimplex *= absK;

        // Compute the gradient
        grad.segment(d*i, d) = gradSimplex;

        // Update the energy
        Ih += absK * G;
    }
    // cout << endl;

    // Now add the constraint regularization
    Ih += 0.5*w*w*( (*DXpU).segment(d*(d+1)*zId, d*(d+1)) - z ).squaredNorm();

    grad += -w*w*(*DXpU).segment(d*(d+1)*zId, d*(d+1)) + w*w*z;

    return Ih;
}

template <int D>
AdaptationFunctional<D>::~AdaptationFunctional() {
}

template class AdaptationFunctional<2>;
template class AdaptationFunctional<3>;