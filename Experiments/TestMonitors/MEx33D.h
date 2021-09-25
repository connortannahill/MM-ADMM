#ifndef M_EX_3_3D_H
#define M_EX_3_3D_H

#include <Eigen/Dense>
#include "../../src/MonitorFunction.h"

template <int D>
class MEx33D : public MonitorFunction<D> {
public:
    void operator()(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &M) override {
        const double PI = 3.141592653589793238462643383;
        M = Eigen::Matrix<double,D,D>::Identity(M.rows(), M.cols());
        if (D == 2) {
            M *= sqrt(0.01/(2.0 + cos(8.0*PI*sqrt(pow(x(0) - 0.5, 2) + pow(x(1) - 0.5, 2)))));
        } else {
            M *= pow(0.01/(2.0 + cos(8.0*PI*sqrt(pow(x(0) - 0.5, 2) + pow(x(1) - 0.5, 2) + pow(x(2) - 0.5, 2)))), 1.0/2.0);
        }
    }
};

template class MEx33D<3>;

#endif