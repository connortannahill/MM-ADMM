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
double HuangFunctional<D>::GC(double (&J)[D][D], double detJ, double (&M)[D][D], double (&x)[D]) {
    double d = (double) D;
    double p = 1.5;
    double theta = 1.0/3.0;
    double deterMinv = 0;
    if(D == 2){
        deterMinv = M[0][0]*M[1][1] - M[0][1]*M[1][0];
    } else {
        double tmp1, tmp2, tmp3;
        tmp1 = M[1][1]*M[2][2] - M[1][2]*M[2][1];
        tmp2 = M[1][0]*M[2][2] - M[1][2]*M[2][0];
        tmp3 = M[1][0]*M[2][1] - M[1][1]*M[2][0];
        deterMinv = M[0][0]*tmp1 - M[0][1]*tmp2 + M[0][2]*tmp3;
    }
    double detM = sqrt(1.0 / deterMinv);
    if(D == 2){
        double I[D][D] = {{J[0][0]*M[0][0]+J[0][1]*M[1][0],J[0][0]*M[0][1]+J[0][1]*M[1][1]},
                            {J[1][0]*M[0][0]+J[1][1]*M[1][0],J[1][0]*M[0][1]+J[1][1]*M[1][1]}};
        // I is J*Minv
        double K[D][D] = {0};
        for(int l = 0; l < D; l++){
            for(int q = 0; q < D; q++){
                double xg = 0;
                for(int o = 0; o < D; o++ ){
                    xg += I[l][o]*J[q][o];
                }
                K[l][q] = xg;
            }
        }
        // K is I*J.T
        double trace = K[0][0] + K[1][1];
    return theta * detM * pow(trace, d*p/2.0)
        + (1.0 - 2.0*theta) * pow(d, d*p/2.0) * detM * pow(detJ/detM, p);
    } else {
        // dimension 3
        double I[D][D] = {0};
        double K[D][D] = {0};
        
        for(int l = 0; l < D; l++){
            for(int q = 0; q < D; q++){
                double xg = 0;
                for(int o = 0; o < D; o++ ){
                    xg += J[l][o]*M[o][q];
                }
                I[l][q] = xg;
            }
        }
        // I is J*Minv
        
        for(int l = 0; l < D; l++){
            for(int q = 0; q < D; q++){
                double xg = 0;
                for(int o = 0; o < D; o++ ){
                    xg += I[l][o]*J[q][o];
                }
                K[l][q] = xg;
            }
        }
        // K is I*J.T
        double trace = K[0][0] + K[1][1]+ K[2][2];
    return theta * detM * pow(trace, d*p/2.0)
        + (1.0 - 2.0*theta) * pow(d, d*p/2.0) * detM * pow(detJ/detM, p);

    }
    //3D under construction
    //return theta * detM * pow((J * M * J.transpose()).trace(), d*p/2.0)
    //    + (1.0 - 2.0*theta) * pow(d, d*p/2.0) * detM * pow(detJ/detM, p);
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
void HuangFunctional<D>::dGdJC(double (&J)[D][D], double detJ, double (&M)[D][D], double (&x)[D], double (&out)[D][D]) {
    double d = (double) D;
    double p = 1.5;
    double theta = 1.0/3.0;
    double deterMinv = 0;
    if(D == 2){
        deterMinv = M[0][0]*M[1][1] - M[0][1]*M[1][0];
    } else {
        double tmp1, tmp2, tmp3;
        tmp1 = M[1][1]*M[2][2] - M[1][2]*M[2][1];
        tmp2 = M[1][0]*M[2][2] - M[1][2]*M[2][0];
        tmp3 = M[1][0]*M[2][1] - M[1][1]*M[2][0];
        deterMinv = M[0][0]*tmp1 - M[0][1]*tmp2 + M[0][2]*tmp3;
    }
    double detM = sqrt(1.0 / deterMinv);
    //double detM = sqrt(1.0 / Minv.determinant());
    
    if(D == 2){
        double I[D][D] = {{J[0][0]*M[0][0]+J[0][1]*M[1][0],J[0][0]*M[0][1]+J[0][1]*M[1][1]},
                            {J[1][0]*M[0][0]+J[1][1]*M[1][0],J[1][0]*M[0][1]+J[1][1]*M[1][1]}};
        // I is J*Minv
        double K[D][D] = {0};
        for(int l = 0; l < D; l++){
            for(int q = 0; q < D; q++){
                double xg = 0;
                for(int o = 0; o < D; o++ ){
                    xg += I[l][o]*J[q][o];
                }
                K[l][q] = xg;
            }
        }
        // K is I*J.T
        double trace = K[0][0] + K[1][1];
        double constant_scalar = d*p*theta*detM * pow(trace, d*p/2.0 - 1);
        for(int l = 0; l < D; l++){
            for(int q = 0; q < D; q++){
                double xg = 0;
                for(int o = 0; o < D; o++ ){
                    xg += M[l][o]*J[q][o];
                }
                out[l][q] = constant_scalar*xg;
            }
        }
        // out is calculated
        return;
    } else {
        // dimension 3
        double I[D][D] = {0};
        double K[D][D] = {0};
        
        for(int l = 0; l < D; l++){
            for(int q = 0; q < D; q++){
                double xg = 0;
                for(int o = 0; o < D; o++ ){
                    xg += J[l][o]*M[o][q];
                }
                I[l][q] = xg;
            }
        }
        // I is J*Minv
        
        for(int l = 0; l < D; l++){
            for(int q = 0; q < D; q++){
                double xg = 0;
                for(int o = 0; o < D; o++ ){
                    xg += I[l][o]*J[q][o];
                }
                K[l][q] = xg;
            }
        }
        // K is I*J.T
        double trace = K[0][0] + K[1][1]+ K[2][2];
        double constant_scalar = d*p*theta*detM * pow(trace, d*p/2.0 - 1);
        for(int l = 0; l < D; l++){
            for(int q = 0; q < D; q++){
                double xg = 0;
                for(int o = 0; o < D; o++ ){
                    xg += M[l][o]*J[q][o];
                }
                out[l][q] = constant_scalar*xg;
            }
        }
        return;
    }
    
    
    //out = d*p*theta*detM * pow((J * Minv * J.transpose()).trace(), d*p/2.0 - 1) * Minv * Jt;
    
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
double HuangFunctional<D>::dGddetC(double (&J)[D][D], double detJ, double (&M)[D][D],
                    double (&x)[D]) {
    double d = (double) D;
    double theta = 1.0/2.0;
    double deterMinv = 0;
    if(D == 2){
        deterMinv = M[0][0]*M[1][1] - M[0][1]*M[1][0];
    } else {
        double tmp1, tmp2, tmp3;
        tmp1 = M[1][1]*M[2][2] - M[1][2]*M[2][1];
        tmp2 = M[1][0]*M[2][2] - M[1][2]*M[2][0];
        tmp3 = M[1][0]*M[2][1] - M[1][1]*M[2][0];
        deterMinv = M[0][0]*tmp1 - M[0][1]*tmp2 + M[0][2]*tmp3;
    }
    double detM = sqrt(1.0 / deterMinv);

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
void HuangFunctional<D>::dGdMC(double (&J)[D][D], double detJ,
                        double (&M)[D][D], double (&x)[D],
                        double (&out)[D][D]) {
    double d = (double) D;
    double p = 1.5;
    double theta = 1.0/3.0;
    double deterMinv = 0;
    if(D == 2){
        deterMinv = M[0][0]*M[1][1] - M[0][1]*M[1][0];
    } else {
        double tmp1, tmp2, tmp3;
        tmp1 = M[1][1]*M[2][2] - M[1][2]*M[2][1];
        tmp2 = M[1][0]*M[2][2] - M[1][2]*M[2][0];
        tmp3 = M[1][0]*M[2][1] - M[1][1]*M[2][0];
        deterMinv = M[0][0]*tmp1 - M[0][1]*tmp2 + M[0][2]*tmp3;
    }
    double detM = sqrt(1.0 / deterMinv);

    //Eigen::Matrix<double,D,D> Jt(J.transpose());
    //double tr = (J * Minv * Jt).trace();
    double tr = 0;
    double MJJM[D][D] = {0};
    if(D == 2){
        double I[D][D] = {{J[0][0]*M[0][0]+J[0][1]*M[1][0],J[0][0]*M[0][1]+J[0][1]*M[1][1]},
                            {J[1][0]*M[0][0]+J[1][1]*M[1][0],J[1][0]*M[0][1]+J[1][1]*M[1][1]}};
        // I is J*Minv
        double K[D][D] = {0};
        for(int l = 0; l < D; l++){
            for(int q = 0; q < D; q++){
                double xg = 0;
                double yg = 0;
                for(int o = 0; o < D; o++ ){
                    xg += I[l][o]*J[q][o];
                    yg += I[o][l]*I[o][q];
                }
                K[l][q] = xg;
                MJJM[l][q] = yg;
            }
        }
        // K is I*J.T
        tr = K[0][0] + K[1][1];
    } else {
        // dimension 3
        double I[D][D] = {0};
        double K[D][D] = {0};
        
        for(int l = 0; l < D; l++){
            for(int q = 0; q < D; q++){
                double xg = 0;
                for(int o = 0; o < D; o++ ){
                    xg += J[l][o]*M[o][q];
                }
                I[l][q] = xg;
            }
        }
        // I is J*Minv
        
        for(int l = 0; l < D; l++){
            for(int q = 0; q < D; q++){
                double xg = 0;
                double yg = 0;
                for(int o = 0; o < D; o++ ){
                    xg += I[l][o]*J[q][o];
                    yg += I[o][l]*I[o][q];
                }
                K[l][q] = xg;
                MJJM[l][q] = yg;
            }
        }
        // K is I*J.T
        tr = K[0][0] + K[1][1]+ K[2][2];
    }
    double c1 = -0.5*theta*d*p * detM * pow(tr, d*p/2.0 - 1);
    double c2 = 0.5*theta * detM * pow(tr, d*p/2.0);
    double c3 = ((0.5-theta)*(1.0-p)*pow(d, d*p/2.0)) * pow(detM, 1-p) * pow(detJ, p);
    
    for(int l = 0; l < D; l++){
        for(int q = 0; q < D; q++){
            out[l][q] = c1*MJJM[l][q] + c2*M[l][q] + c3*M[l][q];
        }
    }
    // part 3: final calculation
    return;
    // out = -0.5*theta*d*p * detM * pow(tr, d*p/2.0 - 1) * Minv.transpose() * Jt * J * Minv 
    //     + 0.5*theta * detM * pow(tr, d*p/2.0) * Minv 
    //     + ((0.5-theta)*(1.0-p)*pow(d, d*p/2.0)) * pow(detM, 1-p) * pow(detJ, p) * Minv;

    // out = - Minv * J.transpose() * J * Minv;
}


template <int D>
void HuangFunctional<D>::dGdX(Eigen::Matrix<double,D,D> &J, double detJ,
            Eigen::Matrix<double,D,D> &Minv, Eigen::Vector<double,D> &x, Eigen::Vector<double,D> &out) {
    out.setZero();
}

template <int D>
void HuangFunctional<D>::dGdXC(double (&J)[D][D], double detJ,
                        double (&M)[D][D], double (&x)[D],
                        double (&out)[D]) {
    for(int i = 0; i < D; ++i){
        out[i] = 0.0;
    }
    return;
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
