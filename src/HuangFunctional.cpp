#include "AdaptationFunctional.h"
#include <iostream>
#include "HuangFunctional.h"
#include <Eigen/Dense>
#include "MonitorFunction.h"
#include <math.h>

#define DIM 2

using namespace std;

HuangFunctional::HuangFunctional(int d, bool boundaryNode, Eigen::MatrixXd &Vc,
        Eigen::MatrixXd &Vp, Eigen::MatrixXi &F,
        Eigen::VectorXd &DXpU, MonitorFunction *m, double w, double theta, double p)
            : AdaptationFunctional(d, boundaryNode, Vc, Vp, F, DXpU, m, w) {
    
    this->p = p;
    this->theta = theta;
    this->boundaryNode = boundaryNode;
    this->w = w;
}

double HuangFunctional::G(Eigen::Matrix2d &J, double detJ, Eigen::Matrix2d &M, Eigen::Vector2d &x) {
    // double sqrtDetM = sqrt(M.determinant());
    // return theta * sqrtDetM * pow((J * M.inverse() * J).trace(), (d*p)/2.0) +
    //     (1 - 2.0*theta)*pow(d, d*p/2.0)*sqrtDetM*pow(detJ/sqrtDetM, p);
    return (J * M.inverse() * J.transpose()).trace();
}

void HuangFunctional::dGdJ(Eigen::Matrix2d &J, double detJ, Eigen::Matrix2d &M, Eigen::Vector2d &x, Eigen::Matrix2d &out) {
    // double sqrtDetM = sqrt(M.determinant());
    // out = d*p*theta*sqrtDetM * pow(((J * M.inverse() * J).trace()), d*p/2.0 - 1) * M.inverse() * J;
    out = 2.0*M.inverse()*J.transpose();
}

double HuangFunctional::dGddet(Eigen::Matrix2d &J, double detJ, Eigen::Matrix2d &M, Eigen::Vector2d &x) {
    // return p*(1.0 - 2.0*theta)*pow(d, d*p/2.0)*pow(M.determinant(), (1.0 - p)/2.0)*pow(detJ, p-1);
    return 0.0;
}

void HuangFunctional::dGdM(Eigen::Matrix2d &J, double detJ, Eigen::Matrix2d &M, Eigen::Vector2d &x, Eigen::Matrix2d &out) {
    // out.setZero();
    // double sqrtDetM = sqrt(M.determinant());

    // out = -theta*(d*p/2.0) * sqrtDetM * pow(((J * M.inverse() * J).trace()), d*p/2.0 - 1) *
    //         M.inverse() * J.transpose() * J * M.inverse() + (theta/2.0) * sqrtDetM *
    //         pow(((J * M.inverse() * J).trace()), d*p/2.0) * M.inverse() +
    //         ((1.0-2.0*theta)*(1.0-p)*pow(d, d*p/2.0))/2.0 * sqrtDetM * pow(detJ/sqrtDetM, p) * M.inverse();

    out = - M.inverse() * J.transpose() * J * M.inverse();
}

void HuangFunctional::dGdX(Eigen::Matrix2d &J, double detJ, Eigen::Matrix2d &M, Eigen::Vector2d &x, Eigen::Vector2d &out) {
    out.setZero();
}

int fact(int n){
     return (n==0) || (n==1) ? 1 : n* fact(n-1);
}

HuangFunctional::HuangFunctional(const HuangFunctional &obj) : AdaptationFunctional(obj) {
    this->d = obj.d;
    this->theta = obj.theta;
    this->w = obj.w;
    this->p = obj.p;
}

/**
 * Method for computing the block gradient of a simplex for the objective function
*/
double HuangFunctional::blockGrad(int zId, Eigen::Vector<double, DIM*(DIM+1)> &z,
            Eigen::Vector<double, DIM*(DIM+1)> &xi,
            Eigen::Vector<double, DIM*(DIM+1)> &grad) {
    double Ih = 0.0;;
    double detFJ;
    Eigen::Vector2d gradSimplex;
    gradSimplex.setZero();
    grad.setZero();

    Eigen::Matrix2d E(Eigen::Matrix2d::Constant(0.0));
    Eigen::Matrix2d FJ(Eigen::Matrix2d::Constant(0.0));
    Eigen::Matrix2d Ehat(Eigen::Matrix2d::Constant(0.0));
    Eigen::Matrix2d M(Eigen::Matrix2d::Constant(0.0));
    Eigen::Vector2d xK(Eigen::Vector2d::Constant(0.0));
    Eigen::Matrix2d vLoc(Eigen::Matrix2d::Constant(0.0));
    Eigen::Matrix2d Mtemp(Eigen::Matrix2d::Constant(0.0));
    Eigen::Vector2d xTemp(Eigen::Vector2d::Constant(0.0));
    Eigen::Matrix2d Einv(Eigen::Matrix2d::Constant(0.0));
    Eigen::Matrix2d dGdJ(Eigen::Matrix2d::Constant(0.0));
    Eigen::Matrix2d Mt0(Eigen::Matrix2d::Constant(0.0));
    Eigen::Matrix2d Mtn(Eigen::Matrix2d::Constant(0.0));
    Eigen::Vector2d dGdX(Eigen::Vector2d::Constant(0.0));
    Eigen::Matrix2d dGdM(Eigen::Matrix2d::Constant(0.0));
    Eigen::Vector2d basisComb(Eigen::Vector2d::Constant(0.0));
    Eigen::Vector<int, DIM+1> ids;

    double dFact = (double) 2.0;
    double absK;
    int d = DIM;

    // Build the edge matrices and the centroid
    xK = z.segment(0, d);
    for (int n = 1; n <= d; n++) {
        xK += z.segment(d*n, d);
    }
    xK /= ((double) d + 1.0);

    // Evaluate the monitor function at the centroid
    M = Mtemp/((double) d + 1);
    (*this->M)(xK, M);

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

        // if (!boundaryNode && computeGrad) {
        this->dGdJ(FJ, detFJ, M, xK, dGdJ);
        double dGddet = this->dGddet(FJ, detFJ, M, xK);
        this->dGdX(FJ, detFJ, M, xK, dGdX);
        this->dGdM(FJ, detFJ, M, xK, dGdM);

        // Form the linear combination of the linear basis derivatives
        xTemp = z.segment(i*d, d);
        (*this->M)(xTemp, Mt0);
        
        basisComb.setZero();
        j = 0;
        for (int n = (i+1)%(d+1); n != i; n = (n+1)%(d+1)) {
            xTemp = z.segment(n*d, d);
            (*this->M)(xTemp, Mtn);
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

    // Now add the constraint regularization
    Ih += 0.5*w*w*( (*DXpU).segment(d*(d+1)*zId, d*(d+1)) - z ).squaredNorm();

    grad += -w*w*(*DXpU).segment(d*(d+1)*zId, d*(d+1)) + w*w*z;

    return Ih;
}