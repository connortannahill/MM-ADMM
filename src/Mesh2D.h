#ifndef MESH_2D_H
#define MESH_2D_H

#include "Assembly.h"
#include <unordered_map>
#include <Eigen/Dense>
#include <vector>
#include "AdaptationFunctional.h"
#include "HuangFunctional.h"
#include "MonitorFunction.h"

#define DIM 2

using namespace std;

class Mesh2D : public Assembly {
public:
    Mesh2D(MonitorFunction *M, unordered_map<string, double> params);
    void prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DxpU, Eigen::VectorXd &z) override;
    // void proxSimplex(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DxpU, Eigen::VectorXd &z) override;
    void updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x) override;
    void copyX(Eigen::VectorXd &tar) override;
    void predictX(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar) override;
    void outputTriangles(const char *fname);
    void outputPoints(const char *fname);
    // double U(int nodeId) override;
    ~Mesh2D();
    // double nablaU(int nodeId) override;
    Eigen::MatrixXd *Vc;
    Eigen::MatrixXd *Vp;
    Eigen::MatrixXi *F;
    Eigen::VectorXi *boundaryMask;
    Eigen::VectorXd *DXpU;
    MonitorFunction *M;
    double newtonOptSimplex(int zId, Eigen::Vector<double, DIM*(DIM+1)> &z,
            Eigen::Vector<double, DIM*(DIM+1)> &xi, int nIter);
    void printDiff();
    int a;
    double tau;
protected:

    // TODO: make this polymorphic
    // vector<HuangFunctional> *I_wx;

    HuangFunctional *I_wx;

    int nx, ny;
    double xa, xb, ya, yb;
    double hx, hy;
    double rho;

    // void buildMassMatrix(Eigen::VectorXd &m) override;
    void buildDMatrix() override;
    void buildWMatrix(double w) override;
};

#endif