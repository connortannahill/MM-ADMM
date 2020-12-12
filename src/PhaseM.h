#ifndef PHASE_M_H
#define PHASE_M_H

#include <Eigen/Dense>
#include "MonitorFunction.h"

class PhaseM : public MonitorFunction {
public:
    void operator()(Eigen::Vector2d &x, Eigen::Matrix2d &M) override;
};

#endif