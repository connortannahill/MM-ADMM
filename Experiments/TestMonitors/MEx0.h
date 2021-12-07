#ifndef M_EX_0_H
#define M_EX_0_H

#include <Eigen/Dense>
#include "../../src/MonitorFunction.h"

template <int D>
class MEx0 : public MonitorFunction<D> {
public:
    void operator()(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &M) override {
        M = Eigen::Matrix<double,D,D>::Identity(M.rows(), M.cols());
        // assert(false);
    }
    
};

template class MEx0<2>;
template class MEx0<3>;

#endif