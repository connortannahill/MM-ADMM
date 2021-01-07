#ifndef ADAPTATION_FUNCTIONAL_H
#define ADAPTATION_FUNCTIONAL_H

#include <Eigen/Dense>
#include "MonitorFunction.h"

class AdaptationFunctional {
protected:
    const Eigen::MatrixXd *Vc;
    const Eigen::MatrixXd *Vp;
    const Eigen::VectorXd *DXpU;
    double w;
    MonitorFunction *M;
    int d;

    virtual double G(Eigen::Matrix2d &J, double det, Eigen::Matrix2d &M, Eigen::Vector2d &x)=0;
    virtual void dGdJ(Eigen::Matrix2d &J, double det, Eigen::Matrix2d &M, Eigen::Vector2d &x, Eigen::Matrix2d &out)=0;
    virtual double dGddet(Eigen::Matrix2d &J, double det, Eigen::Matrix2d &M, Eigen::Vector2d &x)=0;
    virtual void dGdM(Eigen::Matrix2d &J, double det, Eigen::Matrix2d &M, Eigen::Vector2d &x, Eigen::Matrix2d &out)=0;
    virtual void dGdX(Eigen::Matrix2d &J, double det, Eigen::Matrix2d &M, Eigen::Vector2d &x, Eigen::Vector2d &out)=0;
public:
    bool boundaryNode;
    AdaptationFunctional(int d, bool boundaryNode, Eigen::MatrixXd &Vc, Eigen::MatrixXd &Vp,
        Eigen::MatrixXi &F, Eigen::VectorXd &DXpU, MonitorFunction *m, double w);
    AdaptationFunctional(const AdaptationFunctional &obj);
    virtual ~AdaptationFunctional();
private:
};

#endif