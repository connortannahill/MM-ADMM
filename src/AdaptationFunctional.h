#ifndef ADAPTATION_FUNCTIONAL_H
#define ADAPTATION_FUNCTIONAL_H

#include <Eigen/Dense>
#include "MonitorFunction.h"
#include "MeshInterpolator.h"
#include <Eigen/StdVector>

template <int D=-1>
class AdaptationFunctional {
protected:
    const Eigen::MatrixXd *Vc;
    const Eigen::MatrixXd *Vp;
    const Eigen::VectorXd *DXpU;
    double w;
    MonitorFunction<D> *M;
    vector<Eigen::Matrix<double, D, D>> *mPre;

    virtual double G(Eigen::Matrix<double,D,D> &J, double det,
                        Eigen::Matrix<double,D,D> &M,
                        Eigen::Vector<double,D> &x)=0;
    virtual void dGdJ(Eigen::Matrix<double,D,D> &J, double det,
                        Eigen::Matrix<double,D,D> &M, Eigen::Vector<double,D> &x,
                        Eigen::Matrix<double,D,D> &out)=0;
    virtual double dGddet(Eigen::Matrix<double,D,D> &J, double det,
                        Eigen::Matrix<double,D,D> &M, Eigen::Vector<double,D> &x)=0;
    virtual void dGdM(Eigen::Matrix<double,D,D> &J, double det,
                        Eigen::Matrix<double,D,D> &M, Eigen::Vector<double,D> &x,
                        Eigen::Matrix<double,D,D> &out)=0;
    virtual void dGdX(Eigen::Matrix<double,D,D> &J, double det,
                        Eigen::Matrix<double,D,D> &M, Eigen::Vector<double,D> &x,
                        Eigen::Vector<double,D> &out)=0;
public:
    bool boundaryNode;
    AdaptationFunctional(Eigen::MatrixXd &Vc, Eigen::MatrixXd &Vp,
        Eigen::MatrixXi &F, Eigen::VectorXd &DXpU, MonitorFunction<D> *m, double w);
    AdaptationFunctional(const AdaptationFunctional &obj);
    double blockGrad(int zId, Eigen::Vector<double, D*(D+1)> &z,
                Eigen::Vector<double, D*(D+1)> &xi,
                Eigen::Vector<double, D*(D+1)> &grad,
                MeshInterpolator<D> &interp);
    double blockGradC(int zId, Eigen::Vector<double, D*(D+1)> &z,
                Eigen::Vector<double, D*(D+1)> &xi,
                Eigen::Vector<double, D*(D+1)> &grad,
                MeshInterpolator<D> &interp);
    virtual ~AdaptationFunctional();
private:
};

#endif