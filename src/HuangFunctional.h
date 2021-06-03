#ifndef HUANG_FUNCTIONAL_H
#define HUANG_FUNCTIONAL_H

#include <Eigen/Dense>
#include "AdaptationFunctional.h"
// #include "Mesh2D.h"

class Mesh2D;

template <int D>
class HuangFunctional : public AdaptationFunctional<D> {
public:
    HuangFunctional(const HuangFunctional &obj);
    HuangFunctional(Eigen::MatrixXd &Vc, Eigen::MatrixXd &Vp,
        Eigen::MatrixXi &F, Eigen::VectorXd &DXpU, MonitorFunction<D> *m, double w, double theta,
        double p);
    HuangFunctional(MonitorFunction<D> *m,
                    double w, double theta, double p);
    ~HuangFunctional() {};
    double blockGrad(int zId, Eigen::Vector<double, D*(D+1)> &x,
                Eigen::Vector<double, D*(D+1)> &xi,
                Eigen::Vector<double, D*(D+1)> &grad);
protected:
    double theta;
    double p;
    double w;
    double G(Eigen::Matrix<double,D,D> &J, double det, Eigen::Matrix<double,D,D> &M,
                Eigen::Vector<double,D> &x) override;
    double GC(double (&J)[D][D], double det, double (&M)[D][D], double (&x)[D]) override;  
    void dGdJ(Eigen::Matrix<double,D,D> &J, double det, Eigen::Matrix<double,D,D> &M,
                Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &out) override;
    double dGddet(Eigen::Matrix<double,D,D> &J, double det, Eigen::Matrix<double,D,D> &M,
                Eigen::Vector<double,D> &x) override;
    void dGdM(Eigen::Matrix<double,D,D> &J, double det, Eigen::Matrix<double,D,D> &M,
                Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &out) override;
    void dGdX(Eigen::Matrix<double,D,D> &J, double det, Eigen::Matrix<double,D,D> &M,
                Eigen::Vector<double,D> &x, Eigen::Vector<double,D> &out) override;
};

#endif