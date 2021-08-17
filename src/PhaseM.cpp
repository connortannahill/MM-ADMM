#include "PhaseM.h"

#include <Eigen/Dense>
#include <iostream>

using namespace std;

template <int D>
void PhaseM<D>::operator()(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &M) {
    double mu_1 = 20;
    double mu_2 = 20;

    Eigen::Vector<double,D> center(Eigen::Vector<double,D>::Constant(0.5));

    M = Eigen::Matrix<double,D,D>::Identity(M.rows(), M.cols());
    M *= (1 + mu_1/( 1 + mu_2*((x - center).squaredNorm()) ));

    return;

    // double lam1 = 1 + (1.0/cosh(50*(x[0] + x[1] - 1.0)*(x[0] + x[1] - 1.0)));
    // // double lam1 = 1 + (1.0/cosh(50*(x[0] + x[1] - 1.0)*(x[0] + x[1] - 1.0)));
    // // // double lam1 = 1 + 10.0*(1.0/cosh(50*(x[0] + x[1] - 1.0)*(x[0] + x[1] - 1.0)));
    // // // double lam2 = lam1;
    // double lam2 = 1.0/lam1;

    // Eigen::Vector<double,D> v;
    // v << (1.0/sqrt(2.0))*1.0 , (1.0/sqrt(2.0))*1.0; 
    // Eigen::Vector<double,D> vOrth;
    // vOrth << (1.0/sqrt(2.0))*1.0 , -(1.0/sqrt(2.0))*1.0; 

    // Eigen::Matrix<double,D,D> m1(lam1*v*v.transpose());

    // Eigen::Matrix<double,D,D> m2(lam2*vOrth*vOrth.transpose());

    // M = m1 + m2;

    // M = Eigen::Matrix<double,D,D>::Identity(M.rows(), M.cols());
    // double r = sqrt(pow(x(0) - 0.7, 2) + pow(x(1) - 0.5, 2));
    // double theta = atan((x(1) - 0.5) / (x(0) - 0.7) );
    // M *= sqrt(2.0/(5.0 * r * pow(cos(theta - 10.0*r*r), 2.0) + 1.0)  + 1.0);

    // const double PI = 3.141592653589793238462643383;
    // M = Eigen::Matrix<double,D,D>::Identity(M.rows(), M.cols());
    // if (D == 2) {
    //     M *= sqrt(0.01/(2.0 + cos(8.0*PI*sqrt(pow(x(0) - 0.5, 2) + pow(x(1) - 0.5, 2)))));
    // } else {
    //     M *= pow(0.01/(2.0 + cos(8.0*PI*sqrt(pow(x(0) - 0.5, 2) + pow(x(1) - 0.5, 2) + pow(x(2) - 0.5, 2)))), 1.0/2.0);
    // }

    // M = Eigen::Matrix<double,D,D>::Identity(M.rows(), M.cols());
    // double R = 4.5;
    // M += (1.0 / (x(0)*x(0) + sqrt(pow(x(1) - 2.0*R, 2.0) + 1e-8)))*Eigen::Matrix<double,D,D>::Identity(M.rows(), M.cols());

    // double r = sqrt(pow(x(0) - 0.7, 2) + pow(x(1) - 0.5, 2) + pow(x(2) - 0.3, 2));
    // double theta = atan((x(1) - 0.5) / (x(0) - 0.7) );
    // M = Eigen::Matrix<double,D,D>::Identity(M.rows(), M.cols());
    // M *= pow(2.0/(5*r*pow(cos(theta - 10.0*r*r), 2.0) + 1.0)  + 1.0, 1.0/3.0);
}

// Typically only need the 2D and 3D cases. 1D needs to be tested and 4D would break me.
template class PhaseM<2>;
template class PhaseM<3>;
