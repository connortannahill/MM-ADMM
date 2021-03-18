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
    // Eigen::Vector<double,D> centroid;
    Eigen::Matrix<double,D,D> monTemp;
    Eigen::Vector<double,D> xTemp;

    for (int vId = 0; vId < X.rows(); vId++) {
        // Evaluate the montitor function at this vertex
        monTemp.setZero();
        // cout << "extracting the point" << endl;
        xTemp = X.row(vId);
        // cout << "FINISHED extracting the point" << endl;
        // cout << "eval the monitor" << endl;
        (*this)(xTemp, monTemp);
        // cout << "FINISHED eval the monitor" << endl;

        // Place the flattened matrix in the Monitor matrix.
        // cout << "assign to matrix" << endl;
        // cout << "size of monitor row " << M(vId, Eigen::all).size() << endl;
        for (int i = 0; i < D*D; i++) {
            M(vId, i) = monTemp(i/D, i%D);
        }
        // cout << "FINISHED assign to matrix" << endl;
    }
}

/**
 * Function to smooth the monitor function. Just calls the mesh interpolator class
 * to avoid issues with keeping track of memory.
*/
// template <int D>
// void MonitorFunction<D>::smoothMonitor(int nIters, Eigen::MatrixXd &X,
//         Eigen::MatrixXi &F, Eigen::MatrixXd &M, MeshInterpolator<D> &interp) {
//     interp.smoothMonitor(nIters, X, F, M);
// }

// explicit instantiation for each dimension of interest
template class MonitorFunction<2>;
template class MonitorFunction<3>;

