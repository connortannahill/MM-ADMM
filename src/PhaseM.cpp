#include "PhaseM.h"

#include <Eigen/Dense>
#include <iostream>

using namespace std;

template <int D>
void PhaseM<D>::operator()(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &M) {
    double mu_1 = 20;
    double mu_2 = 20;

    // Eigen::Vector<double,D> center(Eigen::Vector<double,D>::Constant(x.size(), 0.5));

    // M = Eigen::Matrix<double,D,D>::Identity(M.rows(), M.cols());
    // M *= (1 + mu_1/( 1 + mu_2*((x - center).squaredNorm()) ));

    double lam1 = 1 + (1.0/cosh(50*(x[0] + x[1] - 1.0)*(x[0] + x[1] - 1.0)));
    // double lam1 = 1 + (1.0/cosh(50*(x[0] + x[1] - 1.0)*(x[0] + x[1] - 1.0)));
    // // double lam1 = 1 + 10.0*(1.0/cosh(50*(x[0] + x[1] - 1.0)*(x[0] + x[1] - 1.0)));
    // // double lam2 = lam1;
    double lam2 = 1.0/lam1;

    Eigen::Vector<double,D> v;
    v << (1.0/sqrt(2.0))*1.0 , (1.0/sqrt(2.0))*1.0; 
    Eigen::Vector<double,D> vOrth;
    vOrth << (1.0/sqrt(2.0))*1.0 , -(1.0/sqrt(2.0))*1.0; 

    Eigen::Matrix<double,D,D> m1(lam1*v*v.transpose());

    Eigen::Matrix<double,D,D> m2(lam2*vOrth*vOrth.transpose());

    M = m1 + m2;
}

// Typically only need the 2D and 3D cases. 1D needs to be tested and 4D would break me.
template class PhaseM<2>;
template class PhaseM<3>;
