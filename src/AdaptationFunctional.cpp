#include "AdaptationFunctional.h"
#include <Eigen/Dense>
#include "MonitorFunction.h"
#include "Mesh2D.h"
#include <vector>
#include <iostream>

#define DIM 2

using namespace std;

AdaptationFunctional::AdaptationFunctional(const AdaptationFunctional &obj) {
    Vc = obj.Vc;
    Vp = obj.Vp;
    DXpU = obj.DXpU;

    this->M = obj.M;
    this->d = obj.d;
    this->boundaryNode = obj.boundaryNode;
}

/**
 * Bring the specified id to the front of the vector via a swap
*/
void swap(int id, Eigen::Vector<int, DIM+1> &ids) {
    for (int i = 0; i < ids.size(); i++) {
        if (ids[i] == id) {
            ids[i] = ids[0];
            ids[0] = id;
            break;
        }
    }
}

AdaptationFunctional::AdaptationFunctional(int d, bool boundaryNode,
        Eigen::MatrixXd &Vc, Eigen::MatrixXd &Vp, Eigen::MatrixXi &F,
        Eigen::VectorXd &DXpU, MonitorFunction *m, double w) {
    
    this->w = w;
    this->d = d;
    this->boundaryNode = boundaryNode;

    this->DXpU = &DXpU;

    this->Vc = &Vc;
    this->Vp = &Vp;
    this->DXpU = &DXpU;

    this->M = m;
}

AdaptationFunctional::~AdaptationFunctional() {
}