// #include "MeshInterpolator.h"
// #include <nanoflann.hpp>
// #include <vector>
// #include <set>
// #include <iostream>
// #include <fstream>
// #include <cstdlib> 
// #include "MeshUtils.h"

// using namespace std;
// using namespace utils;

// const int NUM_RESULTS = 5;

// template<int D>
// MeshInterpolator<D>::MeshInterpolator() {
//     // centroids = new Eigen::MatrixXd(1, D);
//     X = new Eigen::MatrixXd(1, D);
//     F = new Eigen::MatrixXi(1, D);
//     x = new vector<double>();
//     y = new vector<double>();
//     z = new vector<double>();
//     mTemp = new Eigen::MatrixXd(1, D*D);
//     monVals = new Eigen::MatrixXd(1, D*D);
//     monGridVals = new Eigen::MatrixXd(1, D*D);
//     monGridTemp = new Eigen::MatrixXd(1, D*D);
//     connectivity = new vector<vector<int>>();
//     vertexSearchTree
//         = new KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, MeshInterpolator<D>>,
//                 MeshInterpolator<D>, D>
//                 (D /*dim*/, *this, KDTreeSingleIndexAdaptorParams(5 /* max leaf */) );
// }

// template <int D>
// void MeshInterpolator<D>::checkStorage(Eigen::MatrixXd &X, Eigen::MatrixXi &F, bool resize)  {
//     if (F.rows() != (this->F)->rows()) {
//         if (resize) {
//             (this->F)->resize(F.rows(), F.cols());
//         } else {
//             cout << "can not modify the size of the f matrix ahead of this call! call updateMesh" << endl;
//             assert(F.size() != (this->F)->size());
//         }
//     }

//     if (X.rows() != (this->X)->rows()) {
//         if (resize) {
//             (this->mTemp)->resize(X.rows(), D*D);
//             (this->monVals)->resize(X.rows(), D*D);
//             (this->X)->resize(X.rows(), X.cols());
//         } else {
//             assert(X.size() != (this->X)->size());
//         }
//     }
// }

// /** 
//  * Update the location of the vertices.
//  * 
//  * TODO: resize the storage vectors in the case that size(F) or size(X) has changed 
// */
// template <int D>
// void MeshInterpolator<D>::updateMesh(Eigen::MatrixXd &X, Eigen::MatrixXi &F) {
//     // Resize the storage if necissary
//     checkStorage(X, F, true);

//     // Assign these meshes and simplices to the local copies
//     *(this->X) = X;
//     *(this->F) = F;

//     // It is not uneccisary and expensive to call updateMesh unless
//     // the boundary of the domain has changed
//     nx = 2*((int)pow(this->X->size(), 1.0/2.0));
//     ny = 2*((int)pow(this->X->size(), 1.0/2.0));

//     if (D == 2) {
//         nz = 1;
//     } else if (D == 3) {
//         nz = 2*((int)pow(this->X->size(), 1.0/2.0));
//     }

//     monGridVals->resize((nx+1)*(ny+1)*(nz+1), D*D);
//     monGridTemp->resize((nx+1)*(ny+1)*(nz+1), D*D);


//     x->clear();
//     x->resize(nx+1);
//     y->clear();
//     y->resize(ny+1);

//     if (D == 3) {
//         z->clear();
//         z->resize(nz+1);
//     }

//     double xMin = INFINITY;
//     double xMax = -INFINITY;

//     double yMin = INFINITY;
//     double yMax = -INFINITY;

//     double zMin = INFINITY;
//     double zMax = -INFINITY;

//     for (int i = 0; i < this->X->rows(); i++) {
//         xMin = ((*this->X)(i, 0) < xMin) ? (*this->X)(i, 0) : xMin;
//         xMax = ((*this->X)(i, 0) > xMax) ? (*this->X)(i, 0) : xMax;

//         yMin = ((*this->X)(i, 1) < yMin) ? (*this->X)(i, 1) : yMin;
//         yMax = ((*this->X)(i, 1) > yMax) ? (*this->X)(i, 1) : yMax;
        
//         if (D == 3) {
//             zMin = ((*this->X)(i, 2) < zMin) ? (*this->X)(i, 2) : zMin;
//             zMax = ((*this->X)(i, 2) > zMax) ? (*this->X)(i, 2) : zMax;
//         }
//     }

//     utils::linspace(xMin, xMax, nx, *x);
//     utils::linspace(yMin, yMax, ny, *y);
//     if (D == 3)
//         utils::linspace(zMin, zMax, nz, *z);

//     // Build the kdtree search index on the physical mesh and map into the grid
//     vertexSearchTree->buildIndex();
// }

// template <int D>
// void MeshInterpolator<D>::computeBarycentricCoordinates(int simplexId, Eigen::Vector<double, D> &pnt,
//             Eigen::Vector<double, D+1> &bCoords) {

//     if (D == 2) {
//         Eigen::Vector<double, D> x0((*X)((*F)(simplexId, 0), Eigen::all));
//         Eigen::Vector<double, D> x1((*X)((*F)(simplexId, 1), Eigen::all));
//         Eigen::Vector<double, D> x2((*X)((*F)(simplexId, 2), Eigen::all));

//         double det = (x1(1) - x2(1))*(x0(0) - x2(0)) + (x2(0) - x1(0))*(x0(1) - x2(1));
//         bCoords(0) = ( (x1(1) - x2(1))*(pnt(0) - x2(0)) + (x2(0) - x1(0))*(pnt(1) - x2(1)) ) / det;
//         bCoords(1) = ( (x2(1) - x0(1))*(pnt(0) - x2(0)) + (x0(0) - x2(0))*(pnt(1) - x2(1)) ) / det;
//         bCoords(2) = 1 - bCoords(0) - bCoords(1);
//     } else if (D == 3) {
//         Eigen::Matrix<double,D,D> T;
//         for (int i = 0; i < D; i++) {
//             T.col(i) = (*X)((*F)(simplexId,i), Eigen::all) - (*X)((*F)(simplexId,D), Eigen::all);
//         }

//         Eigen::Vector<double, D> temp((*X)((*F)(simplexId, D), Eigen::all));
//         Eigen::Matrix<double,D,D> Tinv(T.inverse());

//         temp = Tinv * (pnt - temp);

//         bCoords.segment(0, D) = temp;

//         bCoords(D) = 1.0;
//         for (int i = 0; i < D; i++) {
//             bCoords(D) -= bCoords(i);
//         }
//     }
// }

// template <int D>
// void MeshInterpolator<D>::nearestNeighGridMap() {
//     std::vector<size_t> ret_index(NUM_RESULTS);
//     std::vector<double> out_dist_sqr(NUM_RESULTS);
//     double query_pt[D];

//     Eigen::Vector<double, D> pnt;
//     if (D == 2) {
//         for (int i = 0; i < nx+1; i++) {
//             for (int j = 0; j < ny+1; j++) {
//                 query_pt[0] = x->at(i);
//                 query_pt[1] = y->at(j);

//                 int numFound = vertexSearchTree->knnSearch(&query_pt[0],
//                     1, &ret_index[0], &out_dist_sqr[0]);
//                 assert(numFound >= 1);
                
//                 (*monGridVals)(j*(nx+1)+i, Eigen::all) = (*monVals)(ret_index.at(0), Eigen::all);
//             }
//         }
//     } else {
//         for (int k = 0; k < nz+1; k++) {
//             for (int i = 0; i < nx+1; i++) {
//                 for (int j = 0; j < ny+1; j++) {
//                     query_pt[0] = x->at(i);
//                     query_pt[1] = y->at(j);
//                     query_pt[2] = z->at(k);

//                     int numFound = vertexSearchTree->knnSearch(&query_pt[0],
//                         1, &ret_index[0], &out_dist_sqr[0]);
//                     assert(numFound >= 1);

//                     (*monGridVals)((nx+1)*(ny+1)*k+i*(nx+1)+j, Eigen::all) = (*monVals)(ret_index.at(0), Eigen::all);
//                 }
//             }
//         }
//     }
// }

// template <int D>
// void MeshInterpolator<D>::interpolateMonitor(MonitorFunction<D> &Mon) {
//     this->Mon = &Mon;
//     // Set up relevant information
//     const int NUM_SMOOTH = 5; // TODO: allow this to be set as a input param
//     cout << "evalling at vertices" << endl;
//     // Mon.evaluateAtVertices(*X, *F, *monVals);;
//     cout << "FINISHED evalling at vertices" << endl;
//     // nearestNeighGridMap();
//     cout << "finished nearestNeigh" << endl;
//     // smoothMonitorGrid(NUM_SMOOTH);
// }

// template <int D>
// void MeshInterpolator<D>::biLinearInterpolation(double x, double y, double xMesh[2],
//                                             double yMesh[2], double *coefs) {
//         double norm = (1/((xMesh[1]-xMesh[0])*(yMesh[1] - yMesh[0])));
//         coefs[0] = norm * (xMesh[1] - x)*(yMesh[1] - y);
//         coefs[1] = norm * (x - xMesh[0])*(yMesh[1] - y);
//         coefs[2] = norm * (xMesh[1] - x)*(y - yMesh[0]);
//         coefs[3] = norm * (x - xMesh[0])*(y - yMesh[0]);
// }

// template <int D>
// inline void MeshInterpolator<D>::evalMonitorOnGrid(Eigen::Vector<double, D> &pnt, Eigen::Matrix<double, D, D> &mVal) {
//     // int xInd = utils::findLimInfMeshPoint(pnt(0), *this->x);
//     // int yInd = utils::findLimInfMeshPoint(pnt(1), *this->y);
//     Eigen::Vector<double, D*D> mTemp;
//     Eigen::Vector<double, D*D> mFlat;

//     if (D == 2) {
//         // double xMesh[2] = {this->x->at(xInd), this->x->at(xInd+1)};
//         // double yMesh[2] = {this->y->at(yInd), this->y->at(yInd+1)};
//         // double x = pnt(0);
//         // double y = pnt(1);
//         // double norm = (1/((xMesh[1]-xMesh[0])*(yMesh[1] - yMesh[0])));
//         // double coefs[4] = {norm * (xMesh[1] - x)*(yMesh[1] - y), 
//         //                    norm * (x - xMesh[0])*(yMesh[1] - y),
//         //                    norm * (xMesh[1] - x)*(y - yMesh[0]),
//         //                    norm * (x - xMesh[0])*(y - yMesh[0])};
//         // // utils::biLinearInterpolation(x(0), x(1), xMesh, yMesh, coefs);
//         // // biLinearInterpolation(x(0), x(1), xMesh, yMesh, coefs);

//         // for (int n = 0; n < D*D; n++) {
//         //     mVal(n / D, n % D) = coefs[0]*(*monGridVals)(yInd*(nx+1)+xInd, n)
//         //                         + coefs[1]*(*monGridVals)(yInd*(nx+1)+xInd+1, n)
//         //                         + coefs[2]*(*monGridVals)((yInd+1)*(nx+1)+xInd, n)
//         //                         + coefs[3]*(*monGridVals)((yInd+1)*(nx+1)+xInd+1, n);
//         // }
//         (*Mon)(pnt, mVal);
//     } else {

//         // int zInd = utils::findLimInfMeshPoint(pnt(2), *this->z);
//         // double coefs[8] = {0.0};
//         // double xMesh[2] = {x->at(xInd), x->at(xInd+1)};
//         // double yMesh[2] = {y->at(yInd), y->at(yInd+1)};
//         // double zMesh[2] = {z->at(zInd), z->at(zInd+1)};
//         // utils::triLinearInterpolation(pnt(0), pnt(1), pnt(2), xMesh, yMesh, zMesh, coefs);

//         // // cout << "Printing coefs" << endl;
//         // // for (int i = 0; i < 8; i++) {
//         // //     cout << coefs[i] << ", ";
//         // // }
//         // // cout << endl;

//         // Eigen::Vector<double, D*D> mFlat(Eigen::Vector<double, D*D>::Zero());

//         // mFlat += coefs[0]*(*monGridVals)(zInd*(nx+1)*(ny+1) + yInd*(nx+1) + xInd, Eigen::all);
//         // mFlat += coefs[1]*(*monGridVals)(zInd*(nx+1)*(ny+1) + yInd*(nx+1) + xInd+1, Eigen::all);
//         // mFlat += coefs[2]*(*monGridVals)(zInd*(nx+1)*(ny+1) + (yInd+1)*(nx+1) + xInd, Eigen::all);
//         // mFlat += coefs[3]*(*monGridVals)(zInd*(nx+1)*(ny+1) + (yInd+1)*(nx+1) + xInd+1, Eigen::all);
//         // mFlat += coefs[4]*(*monGridVals)((zInd+1)*(nx+1)*(ny+1) + yInd*(nx+1) + xInd, Eigen::all);
//         // mFlat += coefs[5]*(*monGridVals)((zInd+1)*(nx+1)*(ny+1) + yInd*(nx+1) + xInd+1, Eigen::all);
//         // mFlat += coefs[6]*(*monGridVals)((zInd+1)*(nx+1)*(ny+1) + (yInd+1)*(nx+1) + xInd, Eigen::all);
//         // mFlat += coefs[7]*(*monGridVals)((zInd+1)*(nx+1)*(ny+1) + (yInd+1)*(nx+1) + xInd+1, Eigen::all);

//         // for (int n = 0; n < D*D; n++) {
//         //     mVal(n / D, n % D) = mFlat(n);
//         // }
//         (*Mon)(pnt, mVal);
//     }
// }

// template <int D>
// void MeshInterpolator<D>::outputGridMesh() {
//     // assert(false);
//     std::ofstream outFile;
//     Eigen::Matrix<double, D, D> mVal;
//     Eigen::Vector<double, D*D> mFlat;
//     Eigen::Vector<double, D> pnt;
//     outFile.open("gridMesh.txt");
//     for (int i = 0; i < nx+1; i ++) {
//         for (int j = 0; j < ny+1; j++) {
//             pnt(0) = x->at(i);
//             pnt(1) = y->at(j);
//             evalMonitorOnGrid(pnt, mVal);
//             // (*Mon)(pnt, mVal);
//             outFile << x->at(i) << ", " << y->at(j) << ", " << mVal.determinant() << endl;
            
//         }
//     }
//     outFile.close();
// }

// template <int D>
// void MeshInterpolator<D>::smoothMonitorGrid(int nIters) { 

//     for (int iter = 0; iter < nIters; iter++) {
//         (*monGridTemp) = (*monGridVals);
//         if (D == 2) {
//             for (int i = 1; i < nx; i++) {
//                 for (int j = 1; j < ny; j++) {
//                     (*monGridVals)(j*(nx+1)+i, Eigen::all) = 0.6*(*monGridTemp)(j*(nx+1) + i, Eigen::all) ;
//                     (*monGridVals)(j*(nx+1)+i, Eigen::all) += 0.1*(*monGridTemp)(j*(nx+1) + i+1, Eigen::all);
//                     (*monGridVals)(j*(nx+1)+i, Eigen::all) += 0.1*(*monGridTemp)(j*(nx+1) + i-1, Eigen::all);
//                     (*monGridVals)(j*(nx+1)+i, Eigen::all) += 0.1*(*monGridTemp)((j+1)*(nx+1) + i, Eigen::all);
//                     (*monGridVals)(j*(nx+1)+i, Eigen::all) += 0.1*(*monGridTemp)((j-1)*(nx+1) + i, Eigen::all);;

//                     // Eigen::Matrix<double, D, D> a;
//                     // for (int n = 0; n < D; n++) {
//                     //     for (int m = 0; m < D; m++) {
//                     //         a(n, m) = (*monGridVals)(j*(nx+1)+i, n*D + m);
//                     //     }
//                     // }

//                     // cout << a.determinant() << endl;
//                 }
//             }
//         } else {
//             double h = 0.4/6; // Distribute between nodes
//             for (int k = 1; k < nz; k++) {
//                 for (int i = 1; i < nx; i++) {
//                     for (int j = 1; j < ny; j++) {
//                         (*monGridVals)((nx+1)*(ny+1)*k + j*(nx+1) + i, Eigen::all) = 0.6*(*monGridTemp)((nx+1)*(ny+1)*k + j*(nx+1) + i, Eigen::all) 
//                             + h*(*monGridTemp)((nx+1)*(ny+1)*k + j*(nx+1) + (i+1), Eigen::all)
//                             + h*(*monGridTemp)((nx+1)*(ny+1)*k + j*(nx+1) + (i-1), Eigen::all) 
//                             + h*(*monGridTemp)((nx+1)*(ny+1)*k + (j+1)*(nx+1) + i, Eigen::all)
//                             + h*(*monGridTemp)((nx+1)*(ny+1)*k + (j-1)*(nx+1) + i, Eigen::all)
//                             + h*(*monGridTemp)((nx+1)*(ny+1)*(k+1) + j*(nx+1) + i, Eigen::all)
//                             + h*(*monGridTemp)((nx+1)*(ny+1)*(k-1) + j*(nx+1) + i, Eigen::all);
//                     }
//                 }
//             }
//         }
//     }
//     // assert(false);
// }

// template <int D>
// MeshInterpolator<D>::~MeshInterpolator() {
//     // delete centroids;
//     delete mTemp;
//     delete vertexSearchTree;
//     delete monVals;

//     delete monGridTemp;
//     delete monGridVals;
//     delete X;
//     delete F;
//     delete x;
//     delete y;
//     delete z;
//     // delete ret_index;
//     // delete out_dist_sqr;
// }

// // explicit instantiation for each dimension of interest
// template class MeshInterpolator<2>;
// template class MeshInterpolator<3>;

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
                (D /*dim*/, *this, KDTreeSingleIndexAdaptorParams(5 /* max leaf */) );
}

template <int D>
void MeshInterpolator<D>::checkStorage(Eigen::MatrixXd &X, Eigen::MatrixXi &F, bool resize)  {
    if (F.rows() != (this->F)->rows()) {
        if (resize) {
            (this->F)->resize(F.rows(), F.cols());
        } else {
            cout << "can not modify the size of the f matrix ahead of this call! call updateMesh" << endl;
            assert(F.size() != (this->F)->size());
        }
    }

    if (X.rows() != (this->X)->rows()) {
        if (resize) {
            (this->mTemp)->resize(X.rows(), D*D);
            (this->monVals)->resize(X.rows(), D*D);
            (this->X)->resize(X.rows(), X.cols());
        } else {
            assert(X.size() != (this->X)->size());
        }
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

    // It is not uneccisary and expensive to call updateMesh unless
    // the boundary of the domain has changed
    nx = 2*((int)pow(this->X->size(), 1.0/2.0));
    ny = 2*((int)pow(this->X->size(), 1.0/2.0));

    if (D == 2) {
        nz = 1;
    } else if (D == 3) {
        nz = 2*((int)pow(this->X->size(), 1.0/2.0));
    }

    monGridVals->resize((nx+1)*(ny+1)*(nz+1), D*D);
    monGridTemp->resize((nx+1)*(ny+1)*(nz+1), D*D);


    x->clear();
    x->resize(nx+1);
    y->clear();
    y->resize(ny+1);

    if (D == 3) {
        z->clear();
        z->resize(nz+1);
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
    std::vector<size_t> ret_index(NUM_RESULTS);
    std::vector<double> out_dist_sqr(NUM_RESULTS);
    double query_pt[D];

    Eigen::Vector<double, D> pnt;
    if (D == 2) {
        for (int i = 0; i < nx+1; i++) {
            for (int j = 0; j < ny+1; j++) {
                query_pt[0] = x->at(i);
                query_pt[1] = y->at(j);

                int numFound = vertexSearchTree->knnSearch(&query_pt[0],
                    1, &ret_index[0], &out_dist_sqr[0]);
                assert(numFound >= 1);
                
                (*monGridVals)(j*(nx+1)+i, Eigen::all) = (*monVals)(ret_index.at(0), Eigen::all);
            }
        }
    } else {
        for (int k = 0; k < nz+1; k++) {
            for (int i = 0; i < nx+1; i++) {
                for (int j = 0; j < ny+1; j++) {
                    query_pt[0] = x->at(i);
                    query_pt[1] = y->at(j);
                    query_pt[2] = z->at(k);

                    int numFound = vertexSearchTree->knnSearch(&query_pt[0],
                        1, &ret_index[0], &out_dist_sqr[0]);
                    assert(numFound >= 1);

                    (*monGridVals)((nx+1)*(ny+1)*k+i*(nx+1)+j, Eigen::all) = (*monVals)(ret_index.at(0), Eigen::all);
                }
            }
        }
    }
}

template <int D>
void MeshInterpolator<D>::interpolateMonitor(MonitorFunction<D> &Mon) {
    this->Mon = &Mon;
    // Set up relevant information
    const int NUM_SMOOTH = 15; // TODO: allow this to be set as a input param
    Mon.evaluateAtVertices(*X, *F, *monVals);;
    nearestNeighGridMap();
    smoothMonitorGrid(NUM_SMOOTH);
}

template <int D>
void MeshInterpolator<D>::biLinearInterpolation(double x, double y, double xMesh[2],
                                            double yMesh[2], double *coefs) {
        double norm = (1/((xMesh[1]-xMesh[0])*(yMesh[1] - yMesh[0])));
        coefs[0] = norm * (xMesh[1] - x)*(yMesh[1] - y);
        coefs[1] = norm * (x - xMesh[0])*(yMesh[1] - y);
        coefs[2] = norm * (xMesh[1] - x)*(y - yMesh[0]);
        coefs[3] = norm * (x - xMesh[0])*(y - yMesh[0]);
}

template <int D>
inline void MeshInterpolator<D>::evalMonitorOnGrid(Eigen::Vector<double, D> &pnt, Eigen::Matrix<double, D, D> &mVal) {
    int xInd = utils::findLimInfMeshPoint(pnt(0), *this->x);
    int yInd = utils::findLimInfMeshPoint(pnt(1), *this->y);
    Eigen::Vector<double, D*D> mTemp;
    Eigen::Vector<double, D*D> mFlat;

    if (D == 2) {
        double xMesh[2] = {this->x->at(xInd), this->x->at(xInd+1)};
        double yMesh[2] = {this->y->at(yInd), this->y->at(yInd+1)};
        double x = pnt(0);
        double y = pnt(1);
        double norm = (1/((xMesh[1]-xMesh[0])*(yMesh[1] - yMesh[0])));
        double coefs[4] = {norm * (xMesh[1] - x)*(yMesh[1] - y), 
                           norm * (x - xMesh[0])*(yMesh[1] - y),
                           norm * (xMesh[1] - x)*(y - yMesh[0]),
                           norm * (x - xMesh[0])*(y - yMesh[0])};
        // utils::biLinearInterpolation(x(0), x(1), xMesh, yMesh, coefs);
        // biLinearInterpolation(x(0), x(1), xMesh, yMesh, coefs);

        for (int n = 0; n < D*D; n++) {
            mVal(n / D, n % D) = coefs[0]*(*monGridVals)(yInd*(nx+1)+xInd, n)
                                + coefs[1]*(*monGridVals)(yInd*(nx+1)+xInd+1, n)
                                + coefs[2]*(*monGridVals)((yInd+1)*(nx+1)+xInd, n)
                                + coefs[3]*(*monGridVals)((yInd+1)*(nx+1)+xInd+1, n);
        }
    } else {
        int zInd = utils::findLimInfMeshPoint(pnt(2), *this->z);
        double coefs[8] = {0.0};
        double xMesh[2] = {x->at(xInd), x->at(xInd+1)};
        double yMesh[2] = {y->at(yInd), y->at(yInd+1)};
        double zMesh[2] = {z->at(zInd), z->at(zInd+1)};
        utils::triLinearInterpolation(pnt(0), pnt(1), pnt(2), xMesh, yMesh, zMesh, coefs);

        // cout << "Printing coefs" << endl;
        // for (int i = 0; i < 8; i++) {
        //     cout << coefs[i] << ", ";
        // }
        // cout << endl;

        Eigen::Vector<double, D*D> mFlat(Eigen::Vector<double, D*D>::Zero());

        mFlat += coefs[0]*(*monGridVals)(zInd*(nx+1)*(ny+1) + yInd*(nx+1) + xInd, Eigen::all);
        mFlat += coefs[1]*(*monGridVals)(zInd*(nx+1)*(ny+1) + yInd*(nx+1) + xInd+1, Eigen::all);
        mFlat += coefs[2]*(*monGridVals)(zInd*(nx+1)*(ny+1) + (yInd+1)*(nx+1) + xInd, Eigen::all);
        mFlat += coefs[3]*(*monGridVals)(zInd*(nx+1)*(ny+1) + (yInd+1)*(nx+1) + xInd+1, Eigen::all);
        mFlat += coefs[4]*(*monGridVals)((zInd+1)*(nx+1)*(ny+1) + yInd*(nx+1) + xInd, Eigen::all);
        mFlat += coefs[5]*(*monGridVals)((zInd+1)*(nx+1)*(ny+1) + yInd*(nx+1) + xInd+1, Eigen::all);
        mFlat += coefs[6]*(*monGridVals)((zInd+1)*(nx+1)*(ny+1) + (yInd+1)*(nx+1) + xInd, Eigen::all);
        mFlat += coefs[7]*(*monGridVals)((zInd+1)*(nx+1)*(ny+1) + (yInd+1)*(nx+1) + xInd+1, Eigen::all);

        for (int n = 0; n < D*D; n++) {
            mVal(n / D, n % D) = mFlat(n);
        }
    }
}

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
            pnt(0) = x->at(i);
            pnt(1) = y->at(j);
            evalMonitorOnGrid(pnt, mVal);
            // (*Mon)(pnt, mVal);
            outFile << x->at(i) << ", " << y->at(j) << ", " << mVal.determinant() << endl;
            
        }
    }
    outFile.close();
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

                    // Eigen::Matrix<double, D, D> a;
                    // for (int n = 0; n < D; n++) {
                    //     for (int m = 0; m < D; m++) {
                    //         a(n, m) = (*monGridVals)(j*(nx+1)+i, n*D + m);
                    //     }
                    // }

                    // cout << a.determinant() << endl;
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
    // assert(false);
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