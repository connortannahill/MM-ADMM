#ifndef HUANG_FUNCTIONAL_H
#define HUANG_FUNCTIONAL_H

#include <Eigen/Dense>
#include "AdaptationFunctional.h"
// #include "Mesh2D.h"
#define DIM 2

class Mesh2D;

class HuangFunctional : public AdaptationFunctional {
public:
    HuangFunctional(const HuangFunctional &obj);
    HuangFunctional(int d, bool boundaryNode, Eigen::MatrixXd &Vc, Eigen::MatrixXd &Vp,
        Eigen::MatrixXi &F, Eigen::VectorXd &DXpU, MonitorFunction *m, double w, double theta,
        double p);
    HuangFunctional(int d, bool boundaryNode, MonitorFunction *m,
                    double w, double theta, double p);
    ~HuangFunctional() {};
    double blockGrad(int zId, Eigen::Vector<double, DIM*(DIM+1)> &x,
                Eigen::Vector<double, DIM*(DIM+1)> &xi,
                Eigen::Vector<double, DIM*(DIM+1)> &grad);
protected:
    double theta;
    double p;
    double w;
    double G(Eigen::Matrix2d &J, double det, Eigen::Matrix2d &M, Eigen::Vector2d &x) override;
    void dGdJ(Eigen::Matrix2d &J, double det, Eigen::Matrix2d &M, Eigen::Vector2d &x, Eigen::Matrix2d &out) override;
    double dGddet(Eigen::Matrix2d &J, double det, Eigen::Matrix2d &M, Eigen::Vector2d &x) override;
    void dGdM(Eigen::Matrix2d &J, double det, Eigen::Matrix2d &M, Eigen::Vector2d &x, Eigen::Matrix2d &out) override;
    void dGdX(Eigen::Matrix2d &J, double det, Eigen::Matrix2d &M, Eigen::Vector2d &x, Eigen::Vector2d &out) override;
};

#endif