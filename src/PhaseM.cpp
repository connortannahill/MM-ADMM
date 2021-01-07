#include "PhaseM.h"

#include <Eigen/Dense>
#include <iostream>

using namespace std;

void PhaseM::operator()(Eigen::Vector2d &x, Eigen::Matrix2d &M) {
    double mu_1 = 20;
    double mu_2 = 20;

    Eigen::Vector2d center(Eigen::Vector2d::Constant(x.size(), 0.5));

    M = Eigen::Matrix2d::Identity(M.rows(), M.cols());
    M *= (1 + mu_1/( 1 + mu_2*((x - center).squaredNorm()) ));

    // // double lam1 = 1 + 10.0*(1.0/cosh(50*(x[0] + x[1] - 1.0)*(x[0] + x[1] - 1.0)));
    // double lam1 = 1 + 10.0*(1.0/cosh(50*(x[0] + x[1] - 1.0)*(x[0] + x[1] - 1.0)));
    // // double lam2 = lam1;
    // double lam2 = 1.0/lam1;

    // Eigen::Vector2d v;
    // v << (1.0/sqrt(2.0))*1.0 , (1.0/sqrt(2.0))*1.0; 
    // Eigen::Vector2d vOrth;
    // vOrth << (1.0/sqrt(2.0))*1.0 , -(1.0/sqrt(2.0))*1.0; 

    // Eigen::Matrix2d m1(lam1*v*v.transpose());

    // Eigen::Matrix2d m2(lam2*vOrth*vOrth.transpose());

    // M = m1 + m2;
}