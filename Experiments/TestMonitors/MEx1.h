#ifndef M_EX_1_H
#define M_EX_1_H

#include <Eigen/Dense>
#include "../../src/MonitorFunction.h"

template <int D>
class MEx1 : public MonitorFunction<D> {
public:
    void operator()(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &M) override {
        double mu_1 = 20;
        double mu_2 = 20;

        Eigen::Vector<double,D> center(Eigen::Vector<double,D>::Constant(0.5));

        M = Eigen::Matrix<double,D,D>::Identity(M.rows(), M.cols());
        M *= (1 + mu_1/( 1 + mu_2*((x - center).squaredNorm()) ));
    }
    
};

template class MEx1<2>;
template class MEx1<3>;

#endif