#include "MeshInterpolator.h"
#include <nanoflann.hpp>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <cstdlib> 
#include "MeshUtils.h"

using namespace std;
using namespace utils;

const int NUM_RESULTS = 5;

template<int D>
MeshInterpolator<D>::MeshInterpolator() {
    // centroids = new Eigen::MatrixXd(1, D);
    X = new Eigen::MatrixXd(1, D);
    F = new Eigen::MatrixXi(1, D);
    x = new vector<double>();
    y = new vector<double>();
    z = new vector<double>();
    mTemp = new Eigen::MatrixXd(1, D*D);
    monVals = new Eigen::MatrixXd(1, D*D);
    monGridVals = new Eigen::MatrixXd(1, D*D);
    monGridTemp = new Eigen::MatrixXd(1, D*D);
    connectivity = new vector<vector<int>>();
    vertexSearchTree
        = new KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, MeshInterpolator<D>>,
                MeshInterpolator<D>, D>
                (3 /*dim*/, *this, KDTreeSingleIndexAdaptorParams(5 /* max leaf */) );
}

template <int D>
void MeshInterpolator<D>::checkStorage(Eigen::MatrixXd &X, Eigen::MatrixXi &F, bool resize)  {
    // bool sizeChanged = false;
    if (F.rows() != (this->F)->rows()) {
        // sizeChanged = true;
        if (resize) {
            (this->F)->resize(F.rows(), F.cols());
        } else {
            cout << "can not modify the size of the f matrix ahead of this call! call updateMesh" << endl;
            assert(F.size() != (this->F)->size());
        }
    }

    if (X.rows() != (this->X)->rows()) {
        // sizeChanged=true;
        if (resize) {
            (this->mTemp)->resize(X.rows(), D*D);
            (this->monVals)->resize(X.rows(), D*D);
            (this->X)->resize(X.rows(), X.cols());
        } else {
            assert(X.size() != (this->X)->size());
        }
    }

    // If the size changed, update the connectivity
    // if (sizeChanged) {
    //     connectivity->clear();

    //     for (int i = 0; i < X.rows(); i++) {
    //         vector<int> neighs;
    //         findNeighbourPoints(i, neighs);

    //         connectivity->push_back(neighs);
    //     }
    // }
}

/** 
 * Update the location of the vertices.
 * 
 * TODO: resize the storage vectors in the case that size(F) or size(X) has changed 
*/
template <int D>
void MeshInterpolator<D>::updateMesh(Eigen::MatrixXd &X, Eigen::MatrixXi &F) {
    // cout << "In updateMesh" << endl;
    // Resize the storage if necissary
    checkStorage(X, F, true);

    // Assign these meshes and simplices to the local copies
    *(this->X) = X;
    *(this->F) = F;

    // It is not uneccisary and expensive to call updateMesh unless
    // the boundary of the domain has changed
    // cout << "Setting up the grid" << endl;
    nx = 2*((int)sqrt(this->X->size()));
    ny = 2*((int)sqrt(this->X->size()));

    if (D == 2) {
        nz = 1;
    } else if (D == 3) {
        nz = 2*((int)sqrt(this->X->size()));
    }

    monGridVals->resize((nx+1)*(ny+1)*(nz+1), D*D);
    monGridTemp->resize((nx+1)*(ny+1)*(nz+1), D*D);


    x->clear();
    x->resize(nx+1);
    y->clear();
    y->resize(ny+1);

    if (D == 3) {
        z->clear();
        z->resize(nz);
    }

    double xMin = INFINITY;
    double xMax = -INFINITY;

    double yMin = INFINITY;
    double yMax = -INFINITY;

    double zMin = INFINITY;
    double zMax = -INFINITY;

    for (int i = 0; i < this->X->rows(); i++) {
        xMin = ((*this->X)(i, 0) < xMin) ? (*this->X)(i, 0) : xMin;
        xMax = ((*this->X)(i, 0) > xMax) ? (*this->X)(i, 0) : xMax;

        yMin = ((*this->X)(i, 1) < yMin) ? (*this->X)(i, 1) : yMin;
        yMax = ((*this->X)(i, 1) > yMax) ? (*this->X)(i, 1) : yMax;
        
        if (D == 3) {
            zMin = ((*this->X)(i, 2) < zMin) ? (*this->X)(i, 2) : zMin;
            zMax = ((*this->X)(i, 2) > zMax) ? (*this->X)(i, 2) : zMax;
        }
    }

    utils::linspace(xMin, xMax, nx, *x);
    utils::linspace(yMin, yMax, ny, *y);
    if (D == 3)
        utils::linspace(zMin, zMax, nz, *z);
    
    // cout << "x = " << endl;
    // for (int i = 0; i < nx+1; i++)  {
    //     cout << x->at(i) << ", ";
    // }
    // cout << endl;
    // assert(false);

    // Build the kdtree search index on the physical mesh and map into the grid
    vertexSearchTree->buildIndex();
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
void MeshInterpolator<D>::nearestNeighGridMap() {
    // cout << "in nearestNeighGridMap" << endl;
    std::vector<size_t> ret_index(NUM_RESULTS);
    std::vector<double> out_dist_sqr(NUM_RESULTS);
    double query_pt[D];

    Eigen::Vector<double, D> pnt;
    if (D == 2) {
        for (int i = 0; i < nx+1; i++) {
            for (int j = 0; j < ny+1; j++) {
                // cout << "(i, j) = (" << i << ", " << j << ")" << endl;
                query_pt[0] = x->at(i);
                // cout << "got x" << endl;
                query_pt[1] = y->at(j);
                // cout << "got y" << endl;

                int numFound = vertexSearchTree->knnSearch(&query_pt[0],
                    1, &ret_index[0], &out_dist_sqr[0]);
                assert(numFound >= 1);
                
                (*monGridVals)(j*(nx+1)+i, Eigen::all) = (*monVals)(ret_index.at(0), Eigen::all);
            }
        }
    } else {
        // int numFound = 0;
        for (int k = 0; k < nz+1; k++) {
            for (int i = 0; i < nx+1; i++) {
                for (int j = 0; j < ny+1; j++) {
                    query_pt[0] = x->at(i);
                    query_pt[1] = y->at(j);
                    query_pt[2] = z->at(k);

                    vertexSearchTree->knnSearch(&query_pt[0],
                        NUM_RESULTS, &ret_index[0], &out_dist_sqr[0]);
                    
                    (*monGridVals)((nx+1)*(ny+1)*k+i*(nx+1)+j, Eigen::all) = (*monVals)(ret_index.at(0), Eigen::all);
                    
                    ret_index.clear();
                    out_dist_sqr.clear();
                }
            }
        }
    }
    // cout << "FINISHED in nearestNeighGridMap" << endl;
}

template <int D>
void MeshInterpolator<D>::interpolateMonitor(MonitorFunction<D> &Mon) {
    // Set up relevant information
    // cout << "Inteprolating monitor funtion" << endl;
    const int NUM_SMOOTH = 3; // TODO: allow this to be set as a input param
    Mon.evaluateAtVertices(*X, *F, *monVals);;
    nearestNeighGridMap();
    smoothMonitorGrid(NUM_SMOOTH);
    outputGridMesh();
    // cout << "FINISHED Inteprolating monitor funtion" << endl;
}

const double CHECK_EPS = 1e-10;

template <int D>
void MeshInterpolator<D>::evalMonitorOnGrid(Eigen::Vector<double, D> &x, Eigen::Matrix<double, D, D> &mVal) {
    // cout << "Evalling monitor on gird" << endl;
    int xInd = utils::findLimInfMeshPoint(x(0), *this->x);
    int yInd = utils::findLimInfMeshPoint(x(1), *this->y);
    // cout << "(x, y) = " << x.transpose() << endl;
    // cout << "xSize = " << this->x->size() << endl;
    // cout << "ySize = " << this->y->size() << endl;
    // cout << "xInd = " << xInd << endl;
    // cout << "yInd = " << yInd << endl;
    // assert(xInd < this->x->size());
    Eigen::Vector<double,D*D> mFlat = Eigen::Vector<double, D*D>::Zero();

    if (D == 2) {
        double coefs[4] = {0.0};
        double xMesh[2] = {this->x->at(xInd), this->x->at(xInd+1)};
        double yMesh[2] = {this->y->at(yInd), this->y->at(yInd+1)};
        utils::biLinearInterpolation(x(0), x(1), xMesh, yMesh, coefs);

        mFlat += coefs[0]*(*monGridVals)(yInd*(nx+1) + xInd, Eigen::all);
        mFlat += coefs[1]*(*monGridVals)(yInd*(nx+1) + xInd+1, Eigen::all);
        mFlat += coefs[2]*(*monGridVals)((yInd+1)*(nx+1) + xInd, Eigen::all);
        mFlat += coefs[3]*(*monGridVals)((yInd+1)*(nx+1) + xInd+1, Eigen::all);
        // cout << "coefs = " << coefs[0] << ", " << coefs[1] << ", " << coefs[2] << ", " << coefs[3] << endl;
    } else {
        int zInd = utils::findLimInfMeshPoint(x(2), *this->z);
        double coefs[8] = {0.0};
        utils::triLinearInterpolation(x(0), x(1), x(2), *this->x, *this->y, *this->z, coefs);

        mFlat += coefs[0]*(*monGridVals)(zInd*(nx+1)*(ny+1) + yInd*(nx+1) + xInd, Eigen::all);
        mFlat += coefs[1]*(*monGridVals)(zInd*(nx+1)*(ny+1) + yInd*(nx+1) + xInd+1, Eigen::all);
        mFlat += coefs[2]*(*monGridVals)(zInd*(nx+1)*(ny+1) + (yInd+1)*(nx+1) + xInd, Eigen::all);
        mFlat += coefs[3]*(*monGridVals)(zInd*(nx+1)*(ny+1) + (yInd+1)*(nx+1) + xInd+1, Eigen::all);
        mFlat += coefs[4]*(*monGridVals)((zInd+1)*(nx+1)*(ny+1) + yInd*(nx+1) + xInd, Eigen::all);
        mFlat += coefs[5]*(*monGridVals)((zInd+1)*(nx+1)*(ny+1) + yInd*(nx+1) + xInd+1, Eigen::all);
        mFlat += coefs[6]*(*monGridVals)((zInd+1)*(nx+1)*(ny+1) + (yInd+1)*(nx+1) + xInd, Eigen::all);
        mFlat += coefs[7]*(*monGridVals)((zInd+1)*(nx+1)*(ny+1) + (yInd+1)*(nx+1) + xInd+1, Eigen::all);
    }

    for (int i = 0; i < D*D; i++) {
        mVal(i/D, i%D) = mFlat(i);
    }
    // cout << "FINISHED Evalling monitor on gird" << endl;
}

// template <int D>
// void MeshInterpolator<D>::evalMonitorOnSimplex(int simplexId, Eigen::Vector<double, D> &x, Eigen::Vector<double,D+1> &bCoords,
//         Eigen::Matrix<double,D,D> &mVal) {
//     // Use the mesh interpolator to find the simplex this point lays on, as well
//     // as its barycentric coordaintes,
//     // Eigen::Vector<double, D+1> bCoords;

//     computeBarycentricCoordinates(simplexId, x, bCoords);

//     // If inside this simplex, stop. Else, perform knn search for proper simplex
//     // bool inTriangle = true;
//     // for (int i = 0; i < D+1; i++)
//     //     inTriangle *= bCoords(i) >= -CHECK_EPS;

//     // // cout << "inTriangle " << inTriangle << endl;

//     // int sId;
//     // if (!inTriangle) {
//         // int sId = evalWithKnn(x, bCoords);
//     // } else {
//     //     sId = simplexId;
//     // }

//     // Now, interpolate the monitor function at this point
//     Eigen::Vector<int, D+1> pntIds((*F)(simplexId, Eigen::all));

//     // Accumulate the monitor function in a flattened array
//     Eigen::Vector<double,D*D> mFlat = Eigen::Vector<double, D*D>::Zero();
//     for (int i = 0; i < D+1; i++) {
//         mFlat += bCoords(i)*(*monVals)(pntIds(i), Eigen::all);
//     }

//     // Now place the value in the output array
//     for (int i = 0; i < D*D; i++) {
//         mVal(i/D, i%D) = mFlat(i);
//     }
// }


template <int D>
void MeshInterpolator<D>::outputGridMesh() {
    // assert(false);
    std::ofstream outFile;
    Eigen::Matrix<double, D, D> mVal;
    Eigen::Vector<double, D*D> mFlat;
    Eigen::Vector<double, D> pnt;
    outFile.open("gridMesh.txt");
    for (int i = 0; i < nx+1; i ++) {
        for (int j = 0; j < ny+1; j++) {
            // mFlat = (*monGridVals)((nx+1)*j+i, Eigen::all);
            pnt(0) = x->at(i);
            pnt(1) = y->at(j);
            evalMonitorOnGrid(pnt, mVal);
            // for (int n = 0; n < D*D; n++) {
            //     mVal(n/D, n%D) = mFlat(n);
            // }
            outFile << x->at(i) << ", " << y->at(j) << ", " << mVal.determinant() << endl;
            
        }
    }
    outFile.close();
    // cout << "finsihed outputting mesh" << endl;
    // assert(false);
}

template <int D>
void MeshInterpolator<D>::smoothMonitorGrid(int nIters) { 

    for (int iter = 0; iter < nIters; iter++) {
        (*monGridTemp) = (*monGridVals);
        if (D == 2) {
            for (int i = 1; i < nx; i++) {
                for (int j = 1; j < ny; j++) {
                    (*monGridVals)(j*(nx+1)+i, Eigen::all) = 0.6*(*monGridTemp)(j*(nx+1) + i, Eigen::all) ;
                    (*monGridVals)(j*(nx+1)+i, Eigen::all) += 0.1*(*monGridTemp)(j*(nx+1) + i+1, Eigen::all);
                    (*monGridVals)(j*(nx+1)+i, Eigen::all) += 0.1*(*monGridTemp)(j*(nx+1) + i-1, Eigen::all);
                    (*monGridVals)(j*(nx+1)+i, Eigen::all) += 0.1*(*monGridTemp)((j+1)*(nx+1) + i, Eigen::all);
                    (*monGridVals)(j*(nx+1)+i, Eigen::all) += 0.1*(*monGridTemp)((j-1)*(nx+1) + i, Eigen::all);;
                }
            }
        } else {
            double h = 0.4/6.0;
            for (int k = 1; k < nz; k++) {
                for (int i = 1; i < nx; i++) {
                    for (int j = 1; j < ny; j++) {
                        (*monGridVals)((nx+1)*(ny+1)*k + j*(nx+1) + i, Eigen::all) = 0.6*(*monGridTemp)((nx+1)*(ny+1)*k + j*(nx+1) + i, Eigen::all) 
                            + h*(*monGridTemp)((nx+1)*(ny+1)*k + j*(nx+1) + (i+1), Eigen::all) + h*(*monGridTemp)((nx+1)*(ny+1)*k + j*(nx+1) + (i-1), Eigen::all) 
                            + h*(*monGridTemp)((nx+1)*(ny+1)*k + (j+1)*(nx+1) + i, Eigen::all) + h*(*monGridTemp)((nx+1)*(ny+1)*k + (j-1)*(nx+1) + i, Eigen::all)
                            + h*(*monGridTemp)((nx+1)*(ny+1)*(k+1) + j*(nx+1) + i, Eigen::all) + h*(*monGridTemp)((nx+1)*(ny+1)*(k-1) + j*(nx+1) + i, Eigen::all);
                    }
                }
            }
        }
    }
}

template <int D>
MeshInterpolator<D>::~MeshInterpolator() {
    // delete centroids;
    delete mTemp;
    delete vertexSearchTree;
    delete monVals;

    delete monGridTemp;
    delete monGridVals;
    delete X;
    delete F;
    delete x;
    delete y;
    delete z;
    // delete ret_index;
    // delete out_dist_sqr;
}

// explicit instantiation for each dimension of interest
template class MeshInterpolator<2>;
template class MeshInterpolator<3>;
