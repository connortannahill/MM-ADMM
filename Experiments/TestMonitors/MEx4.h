#ifndef M_EX_4_H
#define M_EX_4_H

#include <Eigen/Dense>
#include "../../src/MonitorFunction.h"

template <int D>
class MEx4 : public MonitorFunction<D> {
public:
    void operator()(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &M) override {
        double h = 2.0*sqrt(std::numeric_limits<double>::epsilon());
        double eps = 0.01;
        Eigen::Vector<double, 2> grad;

        grad(0) = ((1.0 / (1.0 + exp((x(0)+h + x(1) - 1)/(2.0*eps)))) - (1.0 / (1.0 + exp((x(0)-h + x(1) - 1)/(2.0*eps)))))/(2.0*h);
        grad(1) = ((1.0 / (1.0 + exp((x(0) + x(1)+h - 1)/(2.0*eps)))) - (1.0 / (1.0 + exp((x(0) + x(1)-h - 1)/(2.0*eps)))))/(2.0*h);

        M = Eigen::Matrix<double,D,D>::Identity(M.rows(), M.cols());

        M *= pow(1 + pow(grad.norm(), 2.0), 1.0/4.0);
        // cout << "x = " << x(0) << ", y = " << x(1) << endl;
        // cout << "mult = " << pow(1 + pow(grad.norm(), 2.0), 1.0/2.0) << " det = " << M.determinant() << endl;
    }
};

template class MEx4<2>;
template class MEx4<3>;

#endif