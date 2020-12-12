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
    Floc = new Eigen::MatrixXi(*obj.Floc);
    DXpU = obj.DXpU;

    this->M = obj.M;
    this->d = obj.d;
    this->boundaryNode = obj.boundaryNode;
    this->nodeId = obj.nodeId;
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

AdaptationFunctional::AdaptationFunctional(int d, int nodeId, bool boundaryNode,
        Eigen::MatrixXd &Vc, Eigen::MatrixXd &Vp, Eigen::MatrixXi &F,
        Eigen::VectorXd &DXpU, MonitorFunction *m, double w) {
    
    this->w = w;
    this->d = d;
    this->nodeId = nodeId;
    this->boundaryNode = boundaryNode;

    this->DXpU = &DXpU;

    this->Vc = &Vc;
    this->Vp = &Vp;
    this->DXpU = &DXpU;

    // Get the offsets for all of the simplices connected to this node
    vector<int> inds;
    for (int j = 0; j < F.rows(); j++) {
        if (F(j, 0) == nodeId || F(j, 1) == nodeId || F(j, 2) == nodeId) {
            inds.push_back(j);
        }
    }

    Eigen::Map<Eigen::VectorXi> convMap(&inds[0], inds.size());
    Eigen::VectorXi rowInds(convMap);

    this->Floc = new Eigen::MatrixXi(F(rowInds, Eigen::all));

    // For each row, rearrange to get the target node in the first element
    Eigen::Vector<int, DIM+1> ids;
    for (int i = 0; i < Floc->rows(); i++) {
        ids = Floc->row(i);
        swap(nodeId, ids);
        Floc->row(i) = ids;
    }

    this->M = m;
}

AdaptationFunctional::~AdaptationFunctional() {
    delete Floc;
}