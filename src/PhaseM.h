#ifndef PHASE_M_H
#define PHASE_M_H

#include <Eigen/Dense>
#include "MonitorFunction.h"

template <int D>
class PhaseM : public MonitorFunction<D> {
public:
    void operator()(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &M) override;
};

#endif