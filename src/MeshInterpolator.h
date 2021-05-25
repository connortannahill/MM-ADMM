#ifndef MESH_INTERPOLATOR_H
#define MESH_INTERPOLATOR_H

#include <Eigen/Dense>
#include <nanoflann.hpp>
#include <vector>
#include <iostream>
#include "MonitorFunction.h"

using namespace std;
using namespace nanoflann;

/** 
 * Helper class for interpolating the mesh
*/

template <int D>
class MeshInterpolator {
public:
    MeshInterpolator();
    ~MeshInterpolator();

    Eigen::MatrixXd *centroids;
    Eigen::MatrixXd *mTemp;
    Eigen::MatrixXd *monVals;
    Eigen::MatrixXd *X;
    Eigen::MatrixXi *F;
    vector<vector<int>> *connectivity;

    void updateMesh(Eigen::MatrixXd &X, Eigen::MatrixXi &F);
    int evalWithKnn(Eigen::Vector<double, D> &x, Eigen::Vector<double,D+1> &bCoords);
    void computeBarycentricCoordinates(int simplexId, Eigen::Vector<double, D> &pnt,
            Eigen::Vector<double, D+1> &bCoords);
    void interpolateMonitor(MonitorFunction<D> &Mon);
    void smoothMonitor(int nIters);
    void findNeighbourSimplices(int simplexId, vector<int> neighIds);
    void findNeighbourPoints(int pntId, vector<int> neighPnts);
    void checkStorage(Eigen::MatrixXd &X, Eigen::MatrixXi &F, bool resize);
    void evalMonitorOnSimplex(int simplexId, Eigen::Vector<double, D> &x, Eigen::Vector<double,D+1> &b,
            Eigen::Matrix<double,D,D> &mVal);

    /**
     * Functions and definitions for interfacing with nanoflann KDTree generation.
     * These functions use nearest neighbour search for the closest centroid
    */
    KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, MeshInterpolator<D>>,
            MeshInterpolator<D>, D> *centroidSearchTree;
    // vector<size_t> *ret_index;
    // vector<double> *out_dist_sqr;

    size_t kdtree_get_point_count() const {
		return centroids->rows();
    }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate value, the
    //  "if/else's" are actually solved at compile time.
    // inline double kdtree_get_pt(const size_t idx, const size_t dim) const;
    double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return (*centroids)(idx, dim);
    }

    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const {
        return false;
    }
};

// template <int D=-1>
// class MeshInterpolator {
// public:
//     MeshInterpolator(int nSimplices);
//     ~MeshInterpolator();

//     Eigen::MatrixXd *centroids;

//     void updateMesh(Eigen::MatrixXd &X, Eigen::MatrixXi &F);
//     void eval(Eigen::Vector<double, D> &x, Eigen::MatrixXd &X,
//         Eigen::MatrixXi &F, Eigen::Vector<double,D> *out);
//     void computeBarycentricCoordinates(int simplexId, Eigen::Vector<double, D> &pnt,
//             Eigen::MatrixXd &X, Eigen::MatrixXi &F, Eigen::Vector<double, D+1> &bCoords);
//     void attempt(Eigen::Vector<double,D> &test);


//     /**
//      * Functions and definitions for interfacing with nanoflann KDTree generation.
//      * These functions use nearest neighbour search for the closest centroid
//     */
//     KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, Mesh<D>, Mesh<D>, D> *centroidSearchTree;
//     vector<size_t> *ret_index;
//     vector<num_t> *out_dist_sqr;

//     inline size_t kdtree_get_point_count() const {
// 		return centroids->size();
//     }

// 	// Returns the dim'th component of the idx'th point in the class:
// 	// Since this is inlined and the "dim" argument is typically an immediate value, the
// 	//  "if/else's" are actually solved at compile time.
// 	// inline double kdtree_get_pt(const size_t idx, const size_t dim) const;
// 	inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
//         return (*centroids)(idx, dim);
// 	}

// 	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
// 	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
// 	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
// 	template <class BBOX>
// 	bool kdtree_get_bbox(BBOX& /* bb */) const {
// 		return false;
// 	}

// };

#endif