#include "MeshInterpolator.h"
#include <nanoflann.hpp>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <cstdlib> 

using namespace std;

// const int NUM_RESULTS = 1;

template<int D>
MeshInterpolator<D>::MeshInterpolator() {
    centroids = new Eigen::MatrixXd(1, D);
    X = new Eigen::MatrixXd(1, D);
    F = new Eigen::MatrixXi(1, D);
    mTemp = new Eigen::MatrixXd(1, D*D);
    monVals = new Eigen::MatrixXd(1, D*D);
    connectivity = new vector<vector<int>>();
    centroidSearchTree
        = new KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, MeshInterpolator<D>>,
                MeshInterpolator<D>, D>
                (3 /*dim*/, *this, KDTreeSingleIndexAdaptorParams(5 /* max leaf */) );
}

template <int D>
void MeshInterpolator<D>::checkStorage(Eigen::MatrixXd &X, Eigen::MatrixXi &F, bool resize)  {
    bool sizeChanged = false;
    if (F.rows() != centroids->rows()) {
        sizeChanged = true;
        if (resize) {
            centroids->resize(F.rows(), D);
            (this->F)->resize(F.rows(), F.cols());
        } else {
            cout << "can not modify the size of the f matrix ahead of this call! call updateMesh" << endl;
            assert(F.size() != centroids->size());
        }
    }

    if (X.rows() != (this->mTemp)->rows()) {
        sizeChanged=true;
        if (resize) {
            (this->mTemp)->resize(X.rows(), D*D);
            (this->monVals)->resize(X.rows(), D*D);
            (this->X)->resize(X.rows(), X.cols());
        } else {
            cout << "can not modify the size of the X matrix ahead of this call! call updateMesh" << endl;
            assert(X.size() != (this->mTemp)->size());
        }
    }

    // If the size changed, update the connectivity
    if (sizeChanged) {
	/**
        connectivity->clear();

        for (int i = 0; i < X.rows(); i++) {
            vector<int> neighs;
            findNeighbourPoints(i, neighs);

            connectivity->push_back(neighs);
        }
	*/
    }
}

/** 
 * Update the location of the vertices.
 * 
 * TODO: resize the storage vectors in the case that size(F) or size(X) has changed 
*/
template <int D>
void MeshInterpolator<D>::updateMesh(Eigen::MatrixXd &X, Eigen::MatrixXi &F) {
    // Resize the storage if necissary
    checkStorage(X, F, true);

    // Assign these meshes and simplices to the local copies
    *(this->X) = X;
    *(this->F) = F;

    // Compute the new centroids
    // Eigen::Vector<double,D> temp;
    // Eigen::Vector<double,D> x;
    // for (int i = 0; i < F.rows(); i++) {
    //     temp.setZero();
    //     for (int k = 0; k < D+1; k++) {
    //         x = X(F(i,k), Eigen::all);
    //         temp += x;
    //     }
    //     temp /= (D+1);

    //     for (int j = 0; j < D; j++)  {
    //         (*centroids)(i, j) = temp(j);
    //     }
    // }

    // Update the index
    // centroidSearchTree->buildIndex();
}

template <int D>
void MeshInterpolator<D>::computeBarycentricCoordinates(int simplexId, Eigen::Vector<double, D> &pnt,
            Eigen::Vector<double, D+1> &bCoords) {

    if (D == 2) {
        Eigen::Vector<double, D> x0((*X)((*F)(simplexId, 0), Eigen::all));
        Eigen::Vector<double, D> x1((*X)((*F)(simplexId, 1), Eigen::all));
        Eigen::Vector<double, D> x2((*X)((*F)(simplexId, 2), Eigen::all));

        double det = (x1(1) - x2(1))*(x0(0) - x2(0)) + (x2(0) - x1(0))*(x0(1) - x2(1));
        bCoords(0) = ( (x1(1) - x2(1))*(pnt(0) - x2(0)) + (x2(0) - x1(0))*(pnt(1) - x2(1)) ) / det;
        bCoords(1) = ( (x2(1) - x0(1))*(pnt(0) - x2(0)) + (x0(0) - x2(0))*(pnt(1) - x2(1)) ) / det;
        bCoords(2) = 1 - bCoords(0) - bCoords(1);
    } else if (D == 3) {
        Eigen::Matrix<double,D,D> T;
        for (int i = 0; i < D; i++) {
            T.col(i) = (*X)((*F)(simplexId,i), Eigen::all) - (*X)((*F)(simplexId,D), Eigen::all);
        }

        Eigen::Vector<double, D> temp((*X)((*F)(simplexId, D), Eigen::all));
        Eigen::Matrix<double,D,D> Tinv(T.inverse());

        temp = Tinv * (pnt - temp);

        bCoords.segment(0, D) = temp;

        bCoords(D) = 1.0;
        for (int i = 0; i < D; i++) {
            bCoords(D) -= bCoords(i);
        }
    }
}

template <int D>
void MeshInterpolator<D>::interpolateMonitor(MonitorFunction<D> &Mon) {
    // Set up relevant information
    const int NUM_SMOOTH = 0; // TODO: allow this to be set as a input param
    Mon.evaluateAtVertices(*X, *F, *monVals);;
    smoothMonitor(NUM_SMOOTH);
}

template <int D>
void MeshInterpolator<D>::evalMonitorOnSimplex(int simplexId, Eigen::Vector<double,D> &x,
        Eigen::Matrix<double,D,D> &mVal) {
    // Use the mesh interpolator to find the simplex this point lays on, as well
    // as its barycentric coordaintes,
    Eigen::Vector<double, D+1> bCoords;

    computeBarycentricCoordinates(simplexId, x, bCoords);

    // Now, interpolate the monitor function at this point
    Eigen::Vector<int, D+1> pntIds((*F)(simplexId, Eigen::all));

    // Accumulate the monitor function in a flattened array
    Eigen::Vector<double,D*D> mFlat = Eigen::Vector<double, D*D>::Zero();
    mFlat.setZero();
    for (int i = 0; i < D+1; i++) {
        mFlat += bCoords(i)*(*monVals)(pntIds(i), Eigen::all);
    }

    // Now place the value in the output array
    for (int i = 0; i < D*D; i++) {
        mVal(i/D, i%D) = mFlat(i);
    }
}

/**
 * NOTE: unused
 *
*/
template <int D>
int MeshInterpolator<D>::eval(Eigen::Vector<double, D> &x) {
    assert(false);
    return -1;
    // Do a nearest neighbours search for the nearest centroid
    // std::vector<size_t> ret_index(NUM_RESULTS);
    // std::vector<double> out_dist_sqr(NUM_RESULTS);
    // double query_pt[D];
    // for (int i = 0; i < D; i++) {
    //     query_pt[i] = x(i);
    // }

    // int numFound = centroidSearchTree->knnSearch(&query_pt[0], NUM_RESULTS, &ret_index[0], &out_dist_sqr[0]);

    // // Now, compute the Barycentric coordinates of this point. We check to make sure that
    // // the point is in the correct triangle.
    // bool inTriangle;
    // int simplexId;
    // Eigen::Vector<double, D> roid;
    // const double CHECK_EPS = 1e-8;
    // for (int match = 0; match < numFound; match++) {

    //     // Compute the barycentric coordinates of this point in the matched simplex.
    //     simplexId = ret_index.at(match);
    //     roid = (*centroids)(simplexId, Eigen::all);

    //     // Ensure we are in the simplex with this centroid
    //     inTriangle = true;
        
    //     if (inTriangle) {
    //         // Compute difference in Barycentric reconstruction and the input point
    //         // Eigen::Vector<double, D> test(Eigen::Vector<double, D>::Constant(0.0));
    //         // for (int i = 0; i < D+1; i++) {
    //         //     test += bCoords(i) * (*X)((*F)(simplexId, i), Eigen::all);
    //         // }

    //         // cout << "In compute barycentric coordintes " << bCoords.transpose() << endl;

    //         // cout << "Difference in temp vs x " << (test - x).norm() << endl;
    //         // cout << "err = " << (temp - x).norm() << endl;
    //         break;
    //     }
    // }

    // if (!inTriangle) {
    //     cout << "Error in interpolation! Need to handle this case" << endl;
    //     assert(false);
    // }

    // return simplexId;
}

template <int D>
void MeshInterpolator<D>::findNeighbourSimplices(int simplexId, vector<int> neighIds) {
    return;

    Eigen::Vector<int, D> curSimplexIds((*F)(simplexId, Eigen::all));
    Eigen::Vector<int, D> simplexIds(Eigen::Vector<int, D>::Constant(0));

    set<int> comparisonSet;

    neighIds.clear();

    /* Search through all of the simplices and find all of the neighbours  */
    for (int sId = 0; sId < F->rows(); sId++) {
        if (sId == simplexId) {
            continue;
        }
        
        comparisonSet.clear();
        
        // Assign the current possible neighbour
        simplexIds = (*F)(sId, Eigen::all);

        // Find the intersection between the data of the two elements
        for (int i = 0; i < simplexIds.size(); i++)
            comparisonSet.insert(simplexIds(i));

        for (int i = 0; i < curSimplexIds.size(); i++)
            comparisonSet.insert(curSimplexIds(i));

        // If the intersection is non-empty, push back this id into the matches vector
        if (comparisonSet.size() > 0) {
            neighIds.push_back(sId);
        }
    }
}

/**
 * Find all of the points connected to this one.
 * 
 * TODO: this is very inefficient, but can be sped up using an associative data structure
 *       for each point.
*/
template <int D>
void MeshInterpolator<D>::findNeighbourPoints(int pntId, vector<int> neighPnts) {
    Eigen::Vector<int, D+1> simplexIds(Eigen::Vector<int, D+1>::Constant(0));

    std::set<int> pntSet;

    /* Search through all of the simplices and find all of the simplices containing this point.
        For each of these simplices, keep track of the points connected to this point  */
    // cout << "Finding neigh points" << endl;
    for (int sId = 0; sId < F->rows(); sId++) {
        if (pntId == sId)
            continue;
        
        // Assign the current possible neighbour
        simplexIds = (*F)(sId, Eigen::all);

        // Check whether the point is in this simplex.
        bool inCheck = false;
        for (int i = 0; i < D+1; i++) {
            if (simplexIds(i) == pntId) {
                inCheck = true;
                break;
            }
        }

        // If the point is in the simplex, add its neighbours to the set
        int id;
        if (inCheck) {
            for (int i = 0; i < D+1; i++) {
                id = simplexIds(i);
                if (id != pntId) {
                    pntSet.insert(id);
                }
            }
        }
    }

    // Place the set into the output vector
    neighPnts.clear();
    for (auto pnt = pntSet.begin(); pnt != pntSet.end(); ++pnt) {
        neighPnts.push_back(*pnt);
    }
}

/**
 * Smooth the monitor function using N**2log(N) search, where N is the number
 * of simplices.
 * 
 * TODO: improve this algorithm using a better data structure to find the neighbouring
 *       points.
*/
template <int D>
void MeshInterpolator<D>::smoothMonitor(int nIters) {
    vector<int> neighIds;

    (*this->mTemp) = (*this->monVals);

    double weightPnt;
    double weightNeigh;
    for (int iter = 0; iter < nIters; iter++) {
        for (int pId = 0; pId < monVals->rows(); pId++) {
            // neighIds = connectivity->at(pId);


            // Find all of the neighbours to this simplex
            findNeighbourPoints(pId, neighIds);

            // Moving average weights for then current point
            weightPnt = 0.75;
            weightNeigh = (1.0 - weightPnt) / ((double) neighIds.size());

            // Apply the moving average
            (*monVals)(pId, Eigen::all) = weightPnt*(*(this->mTemp))(pId, Eigen::all);
            for (int i = 0; i < neighIds.size(); i++) {
                (*monVals)(pId, Eigen::all) += weightNeigh*(*(this->mTemp))(neighIds.at(i), Eigen::all);
            }
        }
    }
}

template <int D>
MeshInterpolator<D>::~MeshInterpolator() {
    delete centroids;
    delete mTemp;
    delete centroidSearchTree;
    delete monVals;
    delete X;
    delete F;
    // delete ret_index;
    // delete out_dist_sqr;
}

// explicit instantiation for each dimension of interest
template class MeshInterpolator<2>;
template class MeshInterpolator<3>;
