#ifndef MONITOR_FUNCTION_H
#define MONITOR_FUNCTION_H

#include <Eigen/Dense>

template <int D>
class MonitorFunction {
public:
    virtual void operator()(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &M) = 0;
    virtual ~MonitorFunction() {};
};

#endif