#ifndef MONITOR_FUNCTION_H
#define MONITOR_FUNCTION_H

#include <Eigen/Dense>

class MonitorFunction {
public:
    virtual void operator()(Eigen::Vector2d &x, Eigen::Matrix2d &M) = 0;
    virtual ~MonitorFunction() {};
};

#endif