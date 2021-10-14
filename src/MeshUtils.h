#ifndef MESH_UTILS_H
#define MESH_UTILS_H

#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <map>
#include <unordered_map>
#include "Mesh.h"

using namespace std;

namespace utils {
    /**
     * Uniformly partition [xa, xb] into ns subintervals.
     * Result is stored in x.
    */
    inline void linspace(double xa, double xb, int ns, vector<double> &x) {
        x.resize(ns+1);
        for (int i = 0; i < ns+1; i++) {
            x.at(i) = xa + ((double)i)*(xb - xa)/ns;
        }
    }

    inline void computeCentroid2D(int simplexId, Eigen::MatrixXd *X, Eigen::MatrixXi *F, 
            Eigen::Vector<double, 2> &pnt) {

        Eigen::Vector<double, 2> x0((*X)((*F)(simplexId, 0), Eigen::all));
        Eigen::Vector<double, 2> x1((*X)((*F)(simplexId, 1), Eigen::all));
        Eigen::Vector<double, 2> x2((*X)((*F)(simplexId, 2), Eigen::all));

        pnt = 1.0/3.0 * (x0 + x1 + x2);
    }

   /**
     * Perform a bisection-esque search on a 1D array w_mesh to with indices [0, ..., nw-1] for the maximum mesh
     * point w_i for which w_i < w. Return the i that satisfies this.
    */
    inline int findLimInfMeshPoint(double w, vector<double> &w_mesh)
    {
        uint32_t guess =  (int)((w - w_mesh.at(0))/(w_mesh.at(1) - w_mesh.at(0)));
        if (guess < 0) {
            guess = 0;
        } else if (guess > w_mesh.size()-2) {
            guess = w_mesh.size() - 2;
        }
        return guess;
    }

    inline void biLinearInterpolation(double x, double y, double xMesh[2],
                                            double yMesh[2], double *coefs) {
        double norm = (1/((xMesh[1]-xMesh[0])*(yMesh[1] - yMesh[0])));
        coefs[0] = norm * (xMesh[1] - x)*(yMesh[1] - y);
        coefs[1] = norm * (x - xMesh[0])*(yMesh[1] - y);
        coefs[2] = norm * (xMesh[1] - x)*(y - yMesh[0]);
        coefs[3] = norm * (x - xMesh[0])*(y - yMesh[0]);
    }

    inline void triLinearInterpolation(double x, double y, double z, double xMesh[2],
                                        double yMesh[2], double zMesh[2], 
                                        double coefs[8]) {
        double xd = (x - xMesh[0])/(xMesh[1] - xMesh[0]);
        double yd = (y - yMesh[0])/(yMesh[1] - yMesh[0]);
        double zd = (z - zMesh[0])/(zMesh[1] - zMesh[0]);

        coefs[0] = (1 - xd)*(1 - yd)*(1 - zd);
        coefs[1] = xd*(1 - yd)*(1 - zd);
        coefs[2] = (1 - xd)*yd*(1 - zd);
        coefs[3] = xd*yd*(1 - zd);
        coefs[4] = (1 - xd)*(1 - yd)*zd;
        coefs[5] = xd*(1 - yd)*zd;
        coefs[6] = (1 - xd)*yd*zd;
        coefs[7] = xd*yd*zd;
    }

    template <int D>
    inline void generateUniformRectMesh(unordered_map<string,double> params, Eigen::MatrixXd *Vc,
      Eigen::MatrixXi *F, vector<NodeType> *boundaryMask, NodeType bType) {
          
        int nx = (int) params["nx"];
        int ny = (int) params["ny"];
        int nz = (D == 3) ? (int) params["nz"] : 0;

        int xa = params["xa"];
        int xb = params["xb"];
        int ya = params["ya"];
        int yb = params["yb"];

        int za = (D == 3) ? params["za"] : 0;
        int zb = (D == 3) ? params["zb"] : 0;

        double hx = (xb - xa)/((double)nx);
        double hy = (yb - ya)/((double)ny);
        double hz = (D == 3) ? (zb - za)/((double)nz) : 0;

        int off = 0;

        if (D == 2) {
            for (int j = 0; j <= ny; j++) {
                for (int i = 0; i <= nx; i++) {
                    (*Vc)(off, 0) = hx*i;
                    (*Vc)(off, 1) = hy*j;
                    off++;
                }
            }

            // Append the midoints
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    (*Vc)(off, 0) = hx*i + hx/2.0;
                    (*Vc)(off, 1) = hy*j + hy/2.0;

                    off++;
                }
            }

            int stride = (nx+1) * (ny+1);

            off = 0;
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    // Left
                    (*F)(off, 0) = i            + j*(nx+1);
                    (*F)(off, 1) = stride + i   + j*nx;
                    (*F)(off, 2) = i            + (j+1)*(nx+1);
                    off++;

                    // Top
                    (*F)(off, 0) = stride + i   + j*nx;
                    (*F)(off, 1) = i+1          + (j+1)*(nx+1);
                    (*F)(off, 2) = i            + (j+1)*(nx+1);

                    off++;

                    // Right
                    (*F)(off, 0) = stride + i   + j*nx;
                    (*F)(off, 1) = i+1          + (j+1)*(nx+1);
                    (*F)(off, 2) = i+1          + j*(nx+1);

                    off++;

                    // Bot
                    (*F)(off, 0) = i            + j*(nx+1);
                    (*F)(off, 1) = i+1          + j*(nx+1);
                    (*F)(off, 2) = stride + i   + j*nx;

                    off++;
                }
            }

            for (uint64_t i = 0; i < boundaryMask->size(); i++) {
                boundaryMask->at(i) = NodeType::INTERIOR;
            }

            for (int i = 0; i < (nx+1)*(ny+1); i++) {
                int iOff = i % (nx+1);
                int jOff = i / (ny+1);
                bool boundaryPnt = (iOff == 0) || (iOff == nx) || (jOff == 0) || (jOff == ny);

                if (boundaryPnt) {
                    boundaryMask->at(i) = bType;
                } else {
                    boundaryMask->at(i) = NodeType::INTERIOR;
                }


                // Fix the corners
                if ((iOff == 0 && jOff == 0) || (iOff == nx && jOff == 0) 
                        || (iOff == 0 && jOff == ny) || (iOff == nx && jOff == ny)) {
                    boundaryMask->at(i)  = NodeType::BOUNDARY_FIXED;
                }
            }
        } else if (D == 3) {
            for (int k = 0; k <= nz; k++) {
                for (int j = 0; j <= ny; j++) {
                    for (int i = 0; i <= nx; i++) {
                        (*Vc)(off, 0) = hx*i;
                        (*Vc)(off, 1) = hy*j;
                        (*Vc)(off, 2) = hz*k;

                        off++;
                    }
                }
            }

            // Append the midpoints
            for (int k = 0; k < nz; k++) {
                for (int j = 0; j < ny; j++) {
                    for (int i = 0; i < nx; i++) {
                        (*Vc)(off, 0) = hx*i + hx/2.0;
                        (*Vc)(off, 1) = hy*j + hy/2.0;
                        (*Vc)(off, 2) = hz*k + hz/2.0;

                        off++;
                    }
                }
            }

            int stride = (nx+1) * (ny+1) * (nz+1);

            off = 0;
            for (int k = 0; k < nz; k++) {
                for (int j = 0; j < ny; j++) {
                    for (int i = 0; i < nx; i++) {
                        /* Format: xInd + yInd + zInd */
                        int mid = stride + i + j*nx + k*(nx*ny); 

                        // Bot tets
                        (*F)(off, 0) = i            + j*(nx+1)     + k*(nx+1)*(ny+1);
                        (*F)(off, 1) = i+1          + j*(nx+1)     + k*(nx+1)*(ny+1);
                        (*F)(off, 2) = i+1          + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
                        (*F)(off, 3) = mid;

                        off++;

                        (*F)(off, 0) = i            + j*(nx+1)     + k*(nx+1)*(ny+1);
                        (*F)(off, 1) = i            + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
                        (*F)(off, 2) = i+1          + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
                        (*F)(off, 3) = mid;
                        off++;

                        // Top tets
                        (*F)(off, 0) = i            + j*(nx+1)     + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 1) = i+1          + j*(nx+1)     + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 2) = i+1          + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 3) = mid;
                        off++;

                        (*F)(off, 0) = i            + j*(nx+1)     + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 1) = i            + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 2) = i+1          + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 3) = mid;
                        off++;

                        // Left tets
                        (*F)(off, 0) = i            + j*(nx+1)     + k*(nx+1)*(ny+1);
                        (*F)(off, 1) = i            + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
                        (*F)(off, 2) = i            + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 3) = mid;
                        off++;

                        (*F)(off, 0) = i            + j*(nx+1)     + k*(nx+1)*(ny+1);
                        (*F)(off, 1) = i            + j*(nx+1)     + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 2) = i            + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 3) = mid;
                        off++;

                        // Right tets
                        (*F)(off, 0) = i+1          + j*(nx+1)     + k*(nx+1)*(ny+1);
                        (*F)(off, 1) = i+1          + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
                        (*F)(off, 2) = i+1          + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 3) = mid;
                        off++;

                        (*F)(off, 0) = i+1          + j*(nx+1)     + k*(nx+1)*(ny+1);
                        (*F)(off, 1) = i+1          + j*(nx+1)     + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 2) = i+1          + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 3) = mid;
                        off++;
                        
                        // Back tets
                        (*F)(off, 0) = i            + j*(nx+1)     + k*(nx+1)*(ny+1);
                        (*F)(off, 1) = i+1          + (j)*(nx+1)   + k*(nx+1)*(ny+1);
                        (*F)(off, 2) = i            + j*(nx+1)     + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 3) = mid;
                        off++;

                        (*F)(off, 0) = i+1          + j*(nx+1)     + k*(nx+1)*(ny+1);
                        (*F)(off, 1) = i+1          + j*(nx+1)     + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 2) = i            + j*(nx+1)     + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 3) = mid;
                        off++;

                        // Front tets
                        (*F)(off, 0) = i          + (j+1)*(nx+1)     + k*(nx+1)*(ny+1);
                        (*F)(off, 1) = i+1        + (j+1)*(nx+1)     + k*(nx+1)*(ny+1);
                        (*F)(off, 2) = i          + (j+1)*(nx+1)     + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 3) = mid;
                        off++;

                        (*F)(off, 0) = i+1        + (j+1)*(nx+1)     + k*(nx+1)*(ny+1);
                        (*F)(off, 1) = i+1        + (j+1)*(nx+1)     + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 2) = i          + (j+1)*(nx+1)     + (k+1)*(nx+1)*(ny+1);
                        (*F)(off, 3) = mid;
                        off++;
                    }
                }
            }

            for (int i = 0; i < boundaryMask->size(); i++) {
                boundaryMask->at(i) = NodeType::INTERIOR;
            }

            for (int k = 0; k < nz+1; k++) {
                for (int i = 0; i < (nx+1)*(ny+1); i++) {
                    int iOff = i / (nx+1);
                    int jOff = i % (ny+1);
                    bool boundaryPnt = (iOff == 0) || (iOff == nx)
                        || (jOff == 0) || (jOff == ny) || (k == 0)
                        || (k == nz);
                    int off = k*(nx+1)*(ny+1) + i;

                    if (boundaryPnt) {
                        boundaryMask->at(off) = bType;
                    }

                    // Good code
                    bool corner = (iOff == 0 && jOff == 0)
                                || (iOff == nx && jOff == 0)
                                || (iOff == 0 && jOff == ny)
                                || (iOff == nx && jOff == ny)  // end up
                                || (iOff == 0 && k == 0)
                                || (iOff == nx && k == 0)
                                || (iOff == 0 && k == nz)
                                || (iOff == nx && k == nz)  // end vertical
                                || (k == 0 && jOff == 0)
                                || (k == nz && jOff == 0)
                                || (k == 0 && jOff == ny)
                                || (k == nz && jOff == ny); // end horizontal


                    // Fix the corners
                    if (corner) {
                        boundaryMask->at(off) = NodeType::BOUNDARY_FIXED;
                    }
                }
            }
        }
    }

    inline void removeRow(Eigen::MatrixXi &matrix, unsigned int rowToRemove) {
        unsigned int numRows = matrix.rows()-1;
        unsigned int numCols = matrix.cols();

        if( rowToRemove < numRows )
            matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

        matrix.conservativeResize(numRows,numCols);
    }

    inline double phiSdf(double x, double y, std::function<double(double,double)> phiFun) {
        const double h = sqrt(std::numeric_limits<double>::epsilon());
        double grad[2] = {(phiFun(x+h, y) - phiFun(x-h, y))/(2.0*h), 
                       (phiFun(x, y+h) - phiFun(x, y-h))/(2.0*h)
                      };
        double gradNorm = sqrt(grad[0]*grad[0] + grad[1]*grad[1]);

        return phiFun(x, y)/(gradNorm);
    }

    inline double phiSdf(double x, double y, double z, std::function<double(double,double, double)> phiFun) {
        const double h = sqrt(std::numeric_limits<double>::epsilon());
        double grad[3] = {(phiFun(x+h, y, z) - phiFun(x-h, y, z))/(2.0*h), 
                       (phiFun(x, y+h, z) - phiFun(x, y-h, z))/(2.0*h),
                       (phiFun(x, y, z+h) - phiFun(x, y, z-h))/(2.0*h)
                      };
        double gradNorm = sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);

        return phiFun(x, y, z)/(gradNorm);
    }

    inline void interpolateBoundaryLocation(Eigen::Vector<double, 2> &pnt, std::function<double(double, double)> phiFun,
            std::vector<int> &nVals, std::vector<std::tuple<double, double>> &bb) {

        const double h = 2.0*sqrt(std::numeric_limits<double>::epsilon());

        // Compute the inward normal vector
        Eigen::Vector<double, 2> phiGrad;
        phiGrad[0] = (phiSdf(pnt[0]+h, pnt[1], phiFun) - phiSdf(pnt[0]-h, pnt[1], phiFun))/(2.0*h);
        phiGrad[1] = (phiSdf(pnt[0], pnt[1]+h, phiFun) - phiSdf(pnt[0], pnt[1]-h, phiFun))/(2.0*h);

        phiGrad.normalize();

        pnt = pnt - phiSdf(pnt[0], pnt[1], phiFun) * phiGrad;
    }

    inline void interpolateBoundaryLocation(Eigen::Vector<double, 3> &pnt, std::function<double(double, double, double)> phiFun,
            std::vector<int> &nVals, std::vector<std::tuple<double, double>> &bb) {

        const double h = 2.0*sqrt(std::numeric_limits<double>::epsilon());

        // Compute the inward normal vector
        Eigen::Vector<double, 3> phiGrad;
        phiGrad[0] = (phiSdf(pnt[0]+h, pnt[1], pnt[2], phiFun) - phiSdf(pnt[0]-h, pnt[1], pnt[2], phiFun))/(2.0*h);
        phiGrad[1] = (phiSdf(pnt[0], pnt[1]+h, pnt[2], phiFun) - phiSdf(pnt[0], pnt[1]-h, pnt[2], phiFun))/(2.0*h);
        phiGrad[2] = (phiSdf(pnt[0], pnt[1], pnt[2]+h, phiFun) - phiSdf(pnt[0], pnt[1], pnt[2]-h, phiFun))/(2.0*h);

        phiGrad.normalize();

        pnt = pnt - phiSdf(pnt[0], pnt[1], pnt[2], phiFun) * phiGrad;
    }

    inline void meshFromLevelSetFun(std::function<double(double, double)> phiFun, std::vector<int> &nVals,
        std::vector<std::tuple<double, double>> &bb, Eigen::MatrixXd *Vc, Eigen::MatrixXd *Vp, Eigen::MatrixXi *F, 
        std::vector<NodeType> *boundaryMask, NodeType bType) {
        const double EPS = 1e-12;

        // Create the evaluation meshes
        std::vector<double> x;
        std::vector<double> y;

        const int D = 2;

        int nx = nVals.at(0);
        int ny = nVals.at(1);

        linspace(std::get<0>(bb.at(0)), std::get<1>(bb.at(0)), nx, x);
        linspace(std::get<0>(bb.at(1)), std::get<1>(bb.at(1)), ny, y);

        // Generate the unform rectangular mesh for this domain
        std::unordered_map<std::string,double> params;
        params["nx"] = nx;
        params["ny"] = ny;
        params["xa"] = std::get<0>(bb.at(0));
        params["xb"] = std::get<1>(bb.at(0));
        params["ya"] = std::get<0>(bb.at(1));
        params["yb"] = std::get<1>(bb.at(1));

        Vc->resize((nx+1)*(ny+1) + nx*ny, D);
        Vp->resize((nx+1)*(ny+1) + nx*ny, D);
        F->resize(4*nx*ny, D+1);
        boundaryMask->resize(Vc->rows());

        generateUniformRectMesh<D>(params, Vp, F, boundaryMask, bType);

        for (uint32_t i = 0; i < boundaryMask->size(); i++) {
            boundaryMask->at(i) = NodeType::INTERIOR;
        }

        // Find the value of phi at each centroid so we know which to remove
        std::vector<int> idsToBeRemoved;

        Eigen::Vector<double, 2> x0;
        Eigen::Vector<double, 2> x1;
        Eigen::Vector<double, 2> x2;

        for (int sId = 0; sId < F->rows(); sId++) {
            x0 = (*Vp)((*F)(sId, 0), Eigen::all);
            x1 = (*Vp)((*F)(sId, 1), Eigen::all);
            x2 = (*Vp)((*F)(sId, 2), Eigen::all);

            double phiX0 = phiSdf(x0[0], x0[1], phiFun);
            double phiX1 = phiSdf(x1[0], x1[1], phiFun);
            double phiX2 = phiSdf(x2[0], x2[1], phiFun);

            // If all of the sides are outside the mesh domain, remove the simplex entirely
            if (phiX0 > -EPS && phiX1 > -EPS && phiX2 > -EPS) {
                idsToBeRemoved.push_back(sId);
            }
        }

        // Sort the removed Id's into reverse order for easy removal
        std::sort(idsToBeRemoved.begin(), idsToBeRemoved.end(), std::greater<int>());

        for (auto id = idsToBeRemoved.begin(); id != idsToBeRemoved.end(); ++id) {
            removeRow(*F, *id);
        }

        // Now, find all the used points.
        std::set<int> usedPnts;
        for (int i = 0; i < F->rows(); i++) {
            for (int j = 0; j < F->cols(); j++) {
                usedPnts.insert((*F)(i, j));
            }
        }

        // Interpolate the points onto the boundaries
        Eigen::Vector<double, D> xPnt;
        for (auto pnt = usedPnts.begin(); pnt != usedPnts.end(); ++pnt) {
            xPnt = (*Vp)(*pnt, Eigen::all);

            if (phiSdf(xPnt[0], xPnt[1], phiFun) > -EPS) {
                interpolateBoundaryLocation(xPnt, phiFun, nVals, bb);
                boundaryMask->at(*pnt) = bType;
            }

            (*Vp)(*pnt, Eigen::all) = xPnt;
        }

        *Vc = *Vp;
    }

    inline void meshFromLevelSetFun(std::function<double(double, double, double)> phiFun, std::vector<int> &nVals,
        std::vector<std::tuple<double, double>> &bb, Eigen::MatrixXd *Vc, Eigen::MatrixXd *Vp, Eigen::MatrixXi *F, 
        std::vector<NodeType> *boundaryMask, NodeType bType) {
        const double EPS = 1e-12;

        // Create the evaluation meshes
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;

        int nx = nVals.at(0);
        int ny = nVals.at(1);
        int nz = nVals.at(2);

        linspace(std::get<0>(bb.at(0)), std::get<1>(bb.at(0)), nx, x);
        linspace(std::get<0>(bb.at(1)), std::get<1>(bb.at(1)), ny, y);
        linspace(std::get<0>(bb.at(2)), std::get<1>(bb.at(2)), nz, z);

        // Generate the unform rectangular mesh for this domain
        std::unordered_map<std::string,double> params;
        params["nx"] = nx;
        params["ny"] = ny;
        params["nz"] = nz;
        params["xa"] = std::get<0>(bb.at(0));
        params["xb"] = std::get<1>(bb.at(0));
        params["ya"] = std::get<0>(bb.at(1));
        params["yb"] = std::get<1>(bb.at(1));
        params["za"] = std::get<0>(bb.at(2));
        params["zb"] = std::get<1>(bb.at(2));

        Vc->resize((nx+1)*(ny+1)*(nz+1) + nx*ny*nz, 3);
        Vp->resize((nx+1)*(ny+1)*(nz+1) + nx*ny*nz, 3);
        F->resize(12*nx*ny*nz, 3+1);
        boundaryMask->resize(Vc->rows());

        generateUniformRectMesh<3>(params, Vp, F, boundaryMask, bType);

        for (uint32_t i = 0; i < boundaryMask->size(); i++) {
            boundaryMask->at(i) = NodeType::INTERIOR;
        }

        // Find the value of phi at each centroid so we know which to remove
        std::vector<int> idsToBeRemoved;

        Eigen::Vector<double, 3> x0;
        Eigen::Vector<double, 3> x1;
        Eigen::Vector<double, 3> x2;
        Eigen::Vector<double, 3> x3;

        for (int sId = 0; sId < F->rows(); sId++) {
            x0 = (*Vp)((*F)(sId, 0), Eigen::all);
            x1 = (*Vp)((*F)(sId, 1), Eigen::all);
            x2 = (*Vp)((*F)(sId, 2), Eigen::all);
            x3 = (*Vp)((*F)(sId, 3), Eigen::all);

            double phiX0 = phiSdf(x0[0], x0[1], x0[2], phiFun);
            double phiX1 = phiSdf(x1[0], x1[1], x1[2], phiFun);
            double phiX2 = phiSdf(x2[0], x2[1], x2[2], phiFun);
            double phiX3 = phiSdf(x3[0], x3[1], x3[2], phiFun);

            // If all of the sides are outside the mesh domain, remove the simplex entirely
            if (phiX0 > -EPS && phiX1 > -EPS && phiX2 > -EPS && phiX3 > -EPS) {
                idsToBeRemoved.push_back(sId);
            }
        }

        // Sort the removed Id's into reverse order for easy removal
        std::sort(idsToBeRemoved.begin(), idsToBeRemoved.end(), std::greater<int>());

        for (auto id = idsToBeRemoved.begin(); id != idsToBeRemoved.end(); ++id) {
            removeRow(*F, *id);
        }

        // Now, find all the used points.
        std::set<int> usedPnts;
        for (int i = 0; i < F->rows(); i++) {
            for (int j = 0; j < F->cols(); j++) {
                usedPnts.insert((*F)(i, j));
            }
        }

        // Interpolate the points onto the boundaries
        Eigen::Vector<double, 3> xPnt;
        for (auto pnt = usedPnts.begin(); pnt != usedPnts.end(); ++pnt) {
            xPnt = (*Vp)(*pnt, Eigen::all);

            if (phiSdf(xPnt[0], xPnt[1], xPnt[2], phiFun) > -EPS) {
                interpolateBoundaryLocation(xPnt, phiFun, nVals, bb);
                boundaryMask->at(*pnt) = bType;
            }

            (*Vp)(*pnt, Eigen::all) = xPnt;
        }
    }
}

#endif