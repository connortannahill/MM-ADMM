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
    this->N = F.rows();
}

template <int D>
HuangFunctional<D>::HuangFunctional( Eigen::MatrixXd &Vp, Eigen::MatrixXi &F,
        Eigen::VectorXd &DXpU, MonitorFunction<D> *m, double w, double theta, double p)
            : AdaptationFunctional<D>(Vp, F, DXpU, m, w) {
    
    this->p = p;
    this->theta = theta;
    this->w = w;
    this->N = F.rows();
}

template <int D>
double HuangFunctional<D>::G(Eigen::Matrix<double,D,D> &J, double detJ,
            Eigen::Matrix<double,D,D> &Minv, Eigen::Vector<double,D> &x) {
    double d = (double) D;
    double p = 1.5;
    double theta = 1.0/3.0;
    double detM = sqrt(1.0 / Minv.determinant());
    return theta * detM * pow((J * Minv * J.transpose()).trace(), d*p/2.0)
        + (1.0 - 2.0*theta) * pow(d, d*p/2.0) * detM * pow(detJ/detM, p);

    // return (J * Minv * J.transpose()).trace();
}

template <int D>
void HuangFunctional<D>::dGdJ(Eigen::Matrix<double,D,D> &J, double detJ,
            Eigen::Matrix<double,D,D> &Minv, Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &out) {
    double d = (double) D;
    double p = 1.5;
    double theta = 1.0/3.0;
    double detM = sqrt(1.0 / Minv.determinant());
    Eigen::Matrix<double,D,D> Jt(J.transpose());
    out = d*p*theta*detM * pow((J * Minv * Jt).trace(), d*p/2.0 - 1) * Minv * Jt;
    
    // out = 2.0*Minv*J.transpose();
}

template <int D>
double HuangFunctional<D>::dGddet(Eigen::Matrix<double,D,D> &J, double detJ,
            Eigen::Matrix<double,D,D> &Minv, Eigen::Vector<double,D> &x) {
    double d = (double) D;
    double theta = 1.0/2.0;
    double detM = sqrt(1.0 / Minv.determinant());

    return p*(1.0 - 2.0*theta)*pow(d, (d*p)/2.0)*pow(detM, 1.0 - p)*pow(detJ, p-1);

    // return 0.0;
}

template <int D>
void HuangFunctional<D>::dGdM(Eigen::Matrix<double,D,D> &J, double detJ,
            Eigen::Matrix<double,D,D> &Minv, Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &out) {
    double d = (double) D;
    double p = 1.5;
    double theta = 1.0/3.0;
    double detM = sqrt(1.0 / Minv.determinant());
    Eigen::Matrix<double,D,D> Jt(J.transpose());
    double tr = (J * Minv * Jt).trace();
    out = -0.5*theta*d*p * detM * pow(tr, d*p/2.0 - 1) * Minv.transpose() * Jt * J * Minv 
        + 0.5*theta * detM * pow(tr, d*p/2.0) * Minv 
        + ((0.5-theta)*(1.0-p)*pow(d, d*p/2.0)) * pow(detM, 1-p) * pow(detJ, p) * Minv;

    // out = - Minv * J.transpose() * J * Minv;
}

template <int D>
void HuangFunctional<D>::dGdX(Eigen::Matrix<double,D,D> &J, double detJ,
            Eigen::Matrix<double,D,D> &Minv, Eigen::Vector<double,D> &x, Eigen::Vector<double,D> &out) {
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