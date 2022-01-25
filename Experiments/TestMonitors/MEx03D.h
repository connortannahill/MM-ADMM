#ifndef M_EX_0_3D_H
#define M_EX_0_3D_H

#include <Eigen/Dense>
#include "../../src/MonitorFunction.h"

template <int D>
class MEx13D : public MonitorFunction<D> {
public:
    void operator()(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &M) override {
        M = Eigen::Matrix<double,D,D>::Identity(M.rows(), M.cols());
    }
    
};

template class MEx13D<3>;

#endif