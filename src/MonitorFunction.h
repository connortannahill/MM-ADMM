#ifndef MONITOR_FUNCTION_H
#define MONITOR_FUNCTION_H

#include <Eigen/Dense>
#include <vector>

using namespace std;

template <int D>
class MonitorFunction {
protected:
public:
    virtual void operator()(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &M) = 0;
    virtual ~MonitorFunction() {};
    void evaluateAtVertices(Eigen::MatrixXd &X, Eigen::MatrixXi &F, Eigen::MatrixXd &M);
    void evaluateAtPoint(Eigen::MatrixXd &X, Eigen::MatrixXi &F, Eigen::MatrixXd &M);
    // void smoothMonitor(int nIters, Eigen::MatrixXd &X, Eigen::MatrixXi &F,
    //         Eigen::MatrixXd &M);
    // void interpolateMonitor(MeshInterpolator<D> &interp, Eigen::VectorXd &X,
    //     Eigen::MatrixXi &F, Eigen::MatrixXd &M); // TODO: use friend classes to make this uncallable.
};

#endif