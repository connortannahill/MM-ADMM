#ifndef M_EX_2_3D_H
#define M_EX_2_3D_H

#include <Eigen/Dense>
#include "../../src/MonitorFunction.h"

template <int D>
class MEx23D : public MonitorFunction<D> {
public:
    void operator()(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &M) override {
        double lam1 = 1 + (1.0/cosh(50*(x[0] + x[1] - 1.0)*(x[0] + x[1] - 1.0)));
        double lam2 = 1.0/lam1;

        Eigen::Vector<double,D> v;
        v << (1.0/sqrt(2.0))*1.0 , (1.0/sqrt(2.0))*1.0 , (1.0/sqrt(2.0))*1.0;
        Eigen::Vector<double,D> vOrth;
        vOrth << (1.0/sqrt(2.0))*1.0 , -(1.0/sqrt(2.0))*1.0 , -(1.0/sqrt(2.0))*1.0;

        Eigen::Matrix<double,D,D> m1(lam1*v*v.transpose());
        Eigen::Matrix<double,D,D> m2(lam2*vOrth*vOrth.transpose());
        M = m1 + m2;

        cout << M.determinant() << endl;
    }
};

template class MEx23D<3>;

#endif