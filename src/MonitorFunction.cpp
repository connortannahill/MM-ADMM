#include "MonitorFunction.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "MeshInterpolator.h"
#include <set>

using namespace std;

/**
 * Evaluates the monitor function at the nodes of each simplex
 * 
 * The ith row stores the row-major flattened monitor tensor.
*/
template <int D>
void MonitorFunction<D>::evaluateAtVertices(Eigen::MatrixXd &X, Eigen::MatrixXi &F, Eigen::MatrixXd &M) {
    // Loop through the simplices and evaluate the monitor function at the centroid
    Eigen::Matrix<double,D,D> monTemp;
    Eigen::Vector<double,D> xTemp;

    for (int vId = 0; vId < X.rows(); vId++) {
        // Evaluate the montitor function at this vertex
        monTemp.setZero();
        xTemp = X.row(vId);
        (*this)(xTemp, monTemp);

        // Place the flattened matrix in the Monitor matrix.
        for (int i = 0; i < D*D; i++) {
            M(vId, i) = monTemp(i/D, i%D);
        }
    }
}

// explicit instantiation for each dimension of interest
template class MonitorFunction<2>;
template class MonitorFunction<3>;

