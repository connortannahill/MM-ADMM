#include "MeshInterpolator.h"
#include <nanoflann.hpp>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <cstdlib> 

using namespace std;

const int NUM_RESULTS = 1;

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
void MeshInterpolator<D>::outputStuff() {
    double x, y;

    ofstream outFile;
    const char *fname = "centroid.txt";
    outFile.open(fname);

    // for (auto tup = medialAxisPnts->begin(); tup != medialAxisPnts->end(); ++tup) {
    //     x = get<0>(*tup);
    //     y = get<1>(*tup);

    //     outFile << x << ", ";
    //     outFile << y << endl;
    // }

    // cout << "outputting" << endl;
    Eigen::Matrix<double, D, D> mTemp;
    Eigen::Vector<double, D> centroid;
    // cout << "centroid size (" << centroids->rows() << ", " << centroids->cols() << ")" << endl;
    int nVals = 1000;
    // for (int i = 0; i < X->rows(); i++) {
    //     for (int j = 0; j < D; j++) {
    //         outFile << (*X)(i, j) << ", ";
    //     }

    //     Eigen::Vector<double, D> xTemp(X->row(i));
    //     evalMonitorAtPoint(xTemp, mTemp);

    //     // Eigen::Matrix<double, D, D> mExact;
    //     // for (int j = 0; j < D*D; j++) {
    //     //     mExact(j/D, j%D) = (*monVals)(i, j);
    //     // }
    //     // double err = (mTemp - mExact).norm()/(mExact.norm());

    //     outFile << mTemp.determinant() << endl;
    // }



    // }
    // for (int i = 0; i < centroids->rows(); i++) {
    for (int i = 0; i < nVals; i++) {
        // cout << (*centroids)(i, Eigen::all).cols() << endl;
        // cout << centroid.size() << endl;
        for (int j = 0; j < D; j++) {
            // cout << "assigning" << endl;
            centroid(j) = ((double) rand()) / ((double) RAND_MAX); // (*centroids)(i, j);
        }
        // cout << "out" << endl;
        // centroid = (*centroids)(i, Eigen::all);
        for (int j = 0; j < D; j++) {
            outFile << centroid(j) << ", ";
        }

        // cout << "evalling" << endl;
        assert(false);
        evalMonitorAtPoint(centroid, mTemp);
        // cout << "FINISEHD evalling" << endl;


        outFile << mTemp.determinant() << endl;

        // cout << "centroids done" << endl;

        // for (int j = 0; j < d*d; j++) {
        //     mtemp(j/d, j%d) = (*monvals)(i, j);
        // }
        // outFile << mTemp.determinant() << endl;


    }
    // cout << "FINISHED outputting" << endl;

    outFile.close();
}

template <int D>
void MeshInterpolator<D>::computeBarycentricCoordinates(int simplexId, Eigen::Vector<double, D> &pnt,
            Eigen::Vector<double, D+1> &bCoords) {
    // if (D == 3) {
    //     cout << "DO jnot have 3d case finished" << endl;
    //     assert(false);
    // }

    Eigen::Matrix<double,D,D> T;
    for (int i = 0; i < D; i++) {
        T.col(i) = (*X)((*F)(simplexId,i), Eigen::all) - (*X)((*F)(simplexId,D), Eigen::all);
    }

    Eigen::Vector<double, D> temp((*X)((*F)(simplexId, D), Eigen::all));
    Eigen::Matrix<double,D,D> Tinv(T.inverse());

    // temp = T.lu().solve(pnt - temp);
    // // // // cout << "Solving for the barycenters" << endl;
    // // // // Eigen::Vector<double, D> temp(T.lu().solve(temp2));
    temp = Tinv * (pnt - temp);

    // // bCoords.segment(0, D) = ( T.inverse() ) * (pnt - X(F(simplexId,3), Eigen::all));
    // // cout << "computing the bCoords" << endl;
    bCoords.segment(0, D) = temp;

    // if (D == 2) {
    //     Eigen::Vector<double, D> x0((*X)((*F)(simplexId,0), Eigen::all));
    //     Eigen::Vector<double, D> x1((*X)((*F)(simplexId,1), Eigen::all));
    //     Eigen::Vector<double, D> x2((*X)((*F)(simplexId,2), Eigen::all));

    //     double det = (x1(1) - x2(1)) * (x0(0) - x2(0)) + (x2(0) - x1(0)) * (x0(1) - x2(1));
    //     // // cout << "err = " << det - T.determinant() << endl;
    //     bCoords(0) = ( (x1(1) - x2(1)) * (pnt(0) - x2(0)) + (x2(0) - x1(0)) * (pnt(1) - x2(1)) ) / det;
    //     bCoords(1) = ( (x2(1) - x0(1)) * (pnt(0) - x2(0)) + (x0(0) - x2(0)) * (pnt(1) - x2(1)) ) / det;
    //     bCoords(2) = 1.0 - bCoords(1) - bCoords(0);
    // } else if (D == 3) {
    //     Eigen::Vector<double, D> x0((*X)((*F)(simplexId,0), Eigen::all));
    //     Eigen::Vector<double, D> x1((*X)((*F)(simplexId,1), Eigen::all));
    //     Eigen::Vector<double, D> x2((*X)((*F)(simplexId,2), Eigen::all));
    //     Eigen::Vector<double, D> x3((*X)((*F)(simplexId,3), Eigen::all));

    // }


    bCoords(D) = 1.0;
    for (int i = 0; i < D; i++) {
        bCoords(D) -= bCoords(i);
    }
    // cout << "In compute barycentric coordintes " << bCoords.transpose() << endl;


    // Compute difference in Barycentric reconstruction and the input point
    // Eigen::Vector<double, D> test(Eigen::Vector<double, D>::Constant(0.0));
    // for (int i = 0; i < D+1; i++) {
    //     test += bCoords(i) * (*X)((*F)(simplexId, i), Eigen::all);
    // }

    // cout << "Difference in temp vs x " << (test - pnt).norm() << endl;
    // cout << "finisehd the bCoords" << endl;
}

template <int D>
void MeshInterpolator<D>::interpolateMonitor(MonitorFunction<D> &Mon) {
    // Set up relevant information
    const int NUM_SMOOTH = 0; // TODO: allow this to be set as a input param
    // cout << "eval @ vertices" << endl;
    Mon.evaluateAtVertices(*X, *F, *monVals);;
    // cout << "FINISHED eval @ vertices" << endl;
    // cout << "Smoothing the monitor" << endl;
    smoothMonitor(NUM_SMOOTH);
    // cout << "FINISHED Smoothing the monitor" << endl;

    // Checking whether the interpolation is good
    // Eigen::Vector<double, D> xTemp;
    // Eigen::Matrix<double, D, D> Mtemp;
    // Eigen::Matrix<double, D, D> Mtruu;
    // for (int i = 0; i < X->rows(); i++) {
    //     xTemp = (*X)(i, Eigen::all);
    //     evalMonitorAtPoint(xTemp, Mtemp);

    //     for (int j = 0; j < D*D; j++) {
    //         Mtruu(j/D, j%D) = (*monVals)(i, j);
    //     }

    //     cout << "error = " << (Mtemp - Mtruu).norm() << endl;


    // }

    // outputStuff();
    // assert(false);
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
    Eigen::Vector<double,D*D> mFlat;
    mFlat.setZero();
    for (int i = 0; i < D+1; i++) {
        mFlat += bCoords(i)*(*monVals)(pntIds(i), Eigen::all);
    }

    // Now place the value in the output array
    for (int i = 0; i < D*D; i++) {
        mVal(i/D, i%D) = mFlat(i);
    }
}

template <int D>
void MeshInterpolator<D>::evalMonitorAtPoint(Eigen::Vector<double,D> &x, Eigen::Matrix<double,D,D> &mVal) {
    assert(false);
    // Use the mesh interpolator to find the simplex this point lays on, as well
    // as its barycentric coordaintes,
    // Eigen::Vector<double, D+1> bCoords;

    // Extract the simplex at this point
    int simplexId = this->eval(x);//, bCoords);

    // Evaluate on this simplex
    evalMonitorOnSimplex(simplexId, x, mVal);
}

template <int D>
int MeshInterpolator<D>::eval(Eigen::Vector<double, D> &x) {//, Eigen::Vector<double,D+1> &bCoords) {
    // Do a nearest neighbours search for the nearest centroid
    std::vector<size_t> ret_index(NUM_RESULTS);
    std::vector<double> out_dist_sqr(NUM_RESULTS);
    double query_pt[D];
    for (int i = 0; i < D; i++) {
        query_pt[i] = x(i);
    }

    int numFound = centroidSearchTree->knnSearch(&query_pt[0], NUM_RESULTS, &ret_index[0], &out_dist_sqr[0]);
    
    // cout << "found " << numFound << endl;

    // Now, compute the Barycentric coordinates of this point. We check to make sure that
    // the point is in the correct triangle.
    bool inTriangle;
    int simplexId;
    Eigen::Vector<double, D> roid;
    const double CHECK_EPS = 1e-8;
    for (int match = 0; match < numFound; match++) {

        // Compute the barycentric coordinates of this point in the matched simplex.
        simplexId = ret_index.at(match);
        // cout << "Extracting roid" << endl;
        roid = (*centroids)(simplexId, Eigen::all);
        // cout << "FINSIEHD Extracting roid" << endl;
        // cout << "computing barycenters" << endl;
        // computeBarycentricCoordinates(simplexId, x, bCoords);
        // cout << "in eval" << bCoords.transpose() << endl;
        
        // Ensure we are in the simplex with this centroid
        inTriangle = true;
        // for (int i = 0; i < D+1; i++) {
        //     inTriangle &= (bCoords(i) >= -CHECK_EPS);
        // }
        // assert(inTriangle);
        
        if (inTriangle) {
            // Compute difference in Barycentric reconstruction and the input point
            // Eigen::Vector<double, D> test(Eigen::Vector<double, D>::Constant(0.0));
            // for (int i = 0; i < D+1; i++) {
            //     test += bCoords(i) * (*X)((*F)(simplexId, i), Eigen::all);
            // }

            // cout << "In compute barycentric coordintes " << bCoords.transpose() << endl;

            // cout << "Difference in temp vs x " << (test - x).norm() << endl;
            // cout << "err = " << (temp - x).norm() << endl;
            break;
        }
    }

    if (!inTriangle) {
        cout << "Error in interpolation! Need to handle this case" << endl;
        assert(false);
    }

    return simplexId;
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
        // cout << "assigning ids" << endl;
        simplexIds = (*F)(sId, Eigen::all);
        // cout << "FINISHED assigning ids" << endl;

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
        
        // pntSet.clear();
        
        // Assign the current possible neighbour
        // cout << "extracting simplex w/ id = " << sId << endl;
        simplexIds = (*F)(sId, Eigen::all);
        // assert(false);
        // cout << "FINISHED finding first id" << endl;

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
    // assert(false);

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
            // cout << "fouding neigh pnts" << endl;
            findNeighbourPoints(pId, neighIds);
            // cout << "FINISHED fouding neigh pnts" << endl;

            // Moving average weights for then current point
            weightPnt = 0.5;
            weightNeigh = 0.5 / ((double) neighIds.size());

            // Apply the moving average
            // cout << "computing mon" << endl;
            (*monVals)(pId, Eigen::all) = weightPnt*(*(this->mTemp))(pId, Eigen::all);
            for (int i = 0; i < neighIds.size(); i++) {
                (*monVals)(pId, Eigen::all) += weightNeigh*(*(this->mTemp))(neighIds.at(i), Eigen::all);
            }
            // cout << "FINISHED computing mon" << endl;
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
