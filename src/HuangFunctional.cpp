#include "AdaptationFunctional.h"
#include <iostream>
#include "HuangFunctional.h"
#include <Eigen/Dense>
#include "MonitorFunction.h"
#include <math.h>

using namespace std;

template <int D>
HuangFunctional<D>::HuangFunctional(Eigen::MatrixXd &Vc,
        Eigen::MatrixXd &Vp, Eigen::MatrixXi &F,
        Eigen::VectorXd &DXpU, MonitorFunction<D> *m, double w, double theta, double p)
            : AdaptationFunctional<D>(Vc, Vp, F, DXpU, m, w) {
    
    this->p = p;
    this->theta = theta;
    this->w = w;
}

template <int D>
double HuangFunctional<D>::G(Eigen::Matrix<double,D,D> &J, double detJ,
            Eigen::Matrix<double,D,D> &M, Eigen::Vector<double,D> &x) {
    // double d = (double) D;
    // double p = d / 2.0;
    // double theta = 0.5;
    // double sqrtDetM = sqrt(M.determinant());
    // return theta * sqrtDetM * pow((J * M.inverse() * J).trace(), (d*p)/2.0) +
    //     (1 - 2.0*theta)*pow(d, d*p/2.0)*sqrtDetM*pow(detJ/sqrtDetM, p);

    return (J * M.inverse() * J.transpose()).trace();
}

template <int D>
void HuangFunctional<D>::dGdJ(Eigen::Matrix<double,D,D> &J, double detJ,
            Eigen::Matrix<double,D,D> &M, Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &out) {
    // double d = (double) D;
    // double p = d / 2.0;
    // double theta = 0.5;
    // double sqrtDetM = sqrt(M.determinant());
    // out = d*p*theta*sqrtDetM * pow(((J * M.inverse() * J).trace()), d*p/2.0 - 1) * M.inverse() * J;

    out = 2.0*M.inverse()*J.transpose();
}

template <int D>
double HuangFunctional<D>::dGddet(Eigen::Matrix<double,D,D> &J, double detJ,
            Eigen::Matrix<double,D,D> &M, Eigen::Vector<double,D> &x) {
    // double d = (double) D;
    // double p = d / 2.0;
    // double theta = 0.5;
    // return p*(1.0 - 2.0*theta)*pow(d, d*p/2.0)*pow(M.determinant(), (1.0 - p)/2.0)*pow(detJ, p-1);

    return 0.0;
}

template <int D>
void HuangFunctional<D>::dGdM(Eigen::Matrix<double,D,D> &J, double detJ,
            Eigen::Matrix<double,D,D> &M, Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &out) {
    // out.setZero();
    // double d = (double) D;
    // double p = d / 2.0;
    // double theta = 0.5;
    // double sqrtDetM = sqrt(M.determinant());

    // out = -theta*(d*p/2.0) * sqrtDetM * pow(((J * M.inverse() * J).trace()), d*p/2.0 - 1) *
    //         M.inverse() * J.transpose() * J * M.inverse() + (theta/2.0) * sqrtDetM *
    //         pow(((J * M.inverse() * J).trace()), d*p/2.0) * M.inverse() +
    //         ((1.0-2.0*theta)*(1.0-p)*pow(d, d*p/2.0))/2.0 * sqrtDetM * pow(detJ/sqrtDetM, p) * M.inverse();

    out = - M.inverse() * J.transpose() * J * M.inverse();
}

template <int D>
void HuangFunctional<D>::dGdX(Eigen::Matrix<double,D,D> &J, double detJ,
            Eigen::Matrix<double,D,D> &M, Eigen::Vector<double,D> &x, Eigen::Vector<double,D> &out) {
    out.setZero();
}

int fact(int n){
     return (n==0) || (n==1) ? 1 : n* fact(n-1);
}

template <int D>
HuangFunctional<D>::HuangFunctional(const HuangFunctional &obj) : AdaptationFunctional<D>(obj) {
    this->theta = obj.theta;
    this->w = obj.w;
    this->p = obj.p;
}

template class HuangFunctional<2>;
template class HuangFunctional<3>;
