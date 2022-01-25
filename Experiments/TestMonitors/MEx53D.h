#ifndef M_EX_3D_5_H
#define M_EX_3D_5_H

#include <Eigen/Dense>
#include "../../src/MonitorFunction.h"

template <int D>
class MEx5 : public MonitorFunction<D> {
public:
    double u(double x, double y, double z) {
        double r = sqrt(pow(x - 0.7, 2.0) + pow(y - 0.5, 2.0) + pow(z - 0.5, 2));
        double theta = atan((y - 0.5) / (x - 0.7));
        double psi = atan((z - 0.5) / (x - 0.7))

        return 1.0 + 9.0/(1.0 + 100.0*r*r*pow(cos(theta + psi - 20.0*r*r), 2.0));
    }
    void operator()(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &M) override {
        double h = 2.0*sqrt(std::numeric_limits<double>::epsilon());
        Eigen::Vector<double, 2> grad;
        grad(0) = ( u(x(0)+h, x(1), x(2)) - u(x(0)-h, x(1), x(2)) )/(2.0*h);
        grad(1) = ( u(x(0), x(1)+h, x(2)) - u(x(0), x(1)-h, x(2)) )/(2.0*h);
        grad(1) = ( u(x(0), x(1), x(2)+h) - u(x(0), x(1), x(2)-h) )/(2.0*h);

        M = Eigen::Matrix<double,D,D>::Identity(M.rows(), M.cols());

        M *= pow(1 + pow(grad.norm(), 2.0), 1.0/4.0);
    }
};

template class MEx5<2>;
template class MEx5<3>;

#endif