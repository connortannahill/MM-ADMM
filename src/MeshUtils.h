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

    inline void triLinearInterpolation(double x, double y, double z, vector<double> &xMesh,
                                        vector<double> &yMesh, vector<double> &zMesh, 
                                        double coefs[8]) {
        double xd = (x - xMesh.at(0))/(xMesh.at(1) - xMesh.at(0));
        double yd = (y - yMesh.at(0))/(yMesh.at(1) - yMesh.at(0));
        double zd = (z - zMesh.at(0))/(zMesh.at(1) - zMesh.at(0));

        coefs[0] = (1 - xd)*(1 - yd)*(1 - zd);
        coefs[1] = xd*(1 - yd)*(1 - zd);
        coefs[2] = (1 - xd)*yd*(1 - zd);
        coefs[3] = xd*yd*(1 - zd);
        coefs[4] = (1 - xd)*(1 - yd)*zd;
        coefs[5] = xd*(1 - yd)*zd;
        coefs[6] = (1 - xd)*yd*zd;
        coefs[7] = xd*yd*zd;
    }

    inline void generateUniformRectMesh(unordered_map<string,double> params, Eigen::MatrixXd *Vc,
      Eigen::MatrixXi *F, vector<Mesh<2>::NodeType> *boundaryMask, Mesh<2>::NodeType bType) {
          
        int nx = (int) params["nx"];
        int ny = (int) params["ny"];

        int xa = params["xa"];
        int xb = params["xb"];
        int ya = params["ya"];
        int yb = params["yb"];

        double hx = (xb - xa)/((double)nx);
        double hy = (yb - ya)/((double)ny);

        int off = 0;
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
                (*F)(off, 1) = i            + (j+1)*(nx+1);
                (*F)(off, 2) = stride + i   + j*nx;
                off++;

                // Top
                (*F)(off, 0) = i            + (j+1)*(nx+1);
                (*F)(off, 1) = i+1          + (j+1)*(nx+1);
                (*F)(off, 2) = stride + i   + j*nx;

                off++;

                // Right
                (*F)(off, 0) = i+1          + (j+1)*(nx+1);
                (*F)(off, 1) = i+1          + j*(nx+1);
                (*F)(off, 2) = stride + i   + j*nx;

                off++;

                // Bot
                (*F)(off, 0) = i            + j*(nx+1);
                (*F)(off, 1) = i+1          + j*(nx+1);
                (*F)(off, 2) = stride + i   + j*nx;

                off++;
            }
        }

        for (uint64_t i = 0; i < boundaryMask->size(); i++) {
            boundaryMask->at(i) = Mesh<2>::NodeType::INTERIOR;
        }

        for (int i = 0; i < (nx+1)*(ny+1); i++) {
            int iOff = i % (nx+1);
            int jOff = i / (ny+1);
            bool boundaryPnt = (iOff == 0) || (iOff == nx) || (jOff == 0) || (jOff == ny);

            if (boundaryPnt) {
                boundaryMask->at(i) = bType;
            } else {
                boundaryMask->at(i) = Mesh<2>::NodeType::INTERIOR;
            }


            // Fix the corners
            if ((iOff == 0 && jOff == 0) || (iOff == nx && jOff == 0) 
                    || (iOff == 0 && jOff == ny) || (iOff == nx && jOff == ny)) {
                boundaryMask->at(i)  = Mesh<2>::NodeType::BOUNDARY_FIXED;
            }
        }
    }

    inline void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove) {
        unsigned int numRows = matrix.rows()-1;
        unsigned int numCols = matrix.cols();

        if( rowToRemove < numRows )
            matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

        matrix.conservativeResize(numRows,numCols);
    }

    inline void interpolateBoundaryLocation(Eigen::Vector<double, 2> &pnt, std::function<double(double, double)> phiFun,
            std::vector<int> &nVals, std::vector<std::tuple<double, double>> &bb) {

        double xa = std::get<0>(bb.at(0));
        double xb = std::get<1>(bb.at(0));
        double ya = std::get<0>(bb.at(1));
        double yb = std::get<1>(bb.at(1));

        double hx = (xb - xa)/nVals[0];
        double hy = (yb - ya)/nVals[1];

        // Compute the inward normal vector
        Eigen::Vector<double, 2> phiGrad;
        phiGrad[0] = (phiFun(pnt[0]+hx, pnt[1]) - phiFun(pnt[0]-hx, pnt[1]))/(2.0*hx);
        phiGrad[1] = (phiFun(pnt[0], pnt[1]+hy) - phiFun(pnt[0], pnt[1]-hy))/(2.0*hy);
        phiGrad.normalize();
        Eigen::Vector<double, 2> inNorm = - phiGrad;
        double hMove = sqrt(hx*hx + hy*hy);

        // Take the phi values at these points
        Eigen::Vector<double, 2> bndPnt = {pnt[0]+hMove*inNorm[0], pnt[1]+hMove*inNorm[1]};
        double phiOut = phiFun(pnt[0], pnt[1]);
        double phiIn = phiFun(bndPnt[0], bndPnt[1]);

        // Compute the distance to the interface
        Eigen::Vector<double, 2> distVec = {pnt[0] - bndPnt[0], pnt[1] - bndPnt[1]};
        double d = abs(phiOut/(phiIn - phiOut))*distVec.norm();

        // Finally, compute the interface point location
        pnt = bndPnt + d*inNorm;
    }

    inline void meshFromLevelSetFun(std::function<double(double, double)> phiFun, std::vector<int> &nVals,
        std::vector<std::tuple<double, double>> &bb, Eigen::MatrixXd *Vc, Eigen::MatrixXi *F, 
        vector<Mesh<2>::NodeType> *boundaryMask, Mesh<2>::NodeType bType) {

        // const double h = 2.0*sqrt(std::numeric_limits<double>::epsilon());
        
        // Create the evaluation meshes
        vector<double> x;
        vector<double> y;

        const int D = 2;

        int nx = nVals.at(0);
        int ny = nVals.at(1);

        linspace(std::get<0>(bb.at(0)), std::get<1>(bb.at(0)), nx, x);
        linspace(std::get<0>(bb.at(1)), std::get<1>(bb.at(1)), ny, y);

        // Compute the level set function
        // double phi[ny+1][nx+1];
        // for (uint32_t i = 0; i < x.size(); i++) {
        //     for (uint32_t j = 0; j < y.size(); j++) {
        //         phi[j][i] = phiFun(x.at(i), y.at(j));
        //     }
        // }

        // Generate the unform rectangular mesh for this domain
        std::unordered_map<std::string,double> params;
        params["nx"] = nx;
        params["ny"] = ny;
        params["xa"] = std::get<0>(bb.at(0));
        params["xb"] = std::get<1>(bb.at(0));
        params["ya"] = std::get<0>(bb.at(1));
        params["yb"] = std::get<1>(bb.at(1));

        if (Vc != nullptr) {
            delete Vc;
        }
        if (F != nullptr) {
            delete F;
        }

        Vc = new Eigen::MatrixXd((nx+1)*(ny+1) + nx*ny, D);
        F = new Eigen::MatrixXi(4*nx*ny, D+1);

        generateUniformRectMesh(params, Vc, F, boundaryMask, bType);
        for (uint32_t i = 0; i < boundaryMask->size(); i++) {
            boundaryMask->at(i) = Mesh<2>::NodeType::INTERIOR;
        }

        // Find the value of phi at each centroid so we know which to remove
        std::vector<int> idsToBeRemoved;

        Eigen::Vector<double, 2> x0;
        Eigen::Vector<double, 2> x1;
        Eigen::Vector<double, 2> x2;

        for (int sId = 0; sId < F->size(); sId++) {
            x0 = (*Vc)((*F)(sId, 0), Eigen::all);
            x1 = (*Vc)((*F)(sId, 1), Eigen::all);
            x2 = (*Vc)((*F)(sId, 2), Eigen::all);

            double phiX0 = phiFun(x0[0], x0[1]);
            double phiX1 = phiFun(x1[0], x1[1]);
            double phiX2 = phiFun(x2[0], x2[1]);

            // If all of the sides are outside the mesh domain, remove the simplex entirely
            if (phiX0 < 0.0 && phiX1 < 0.0 && phiX2 < 0.0) {
                idsToBeRemoved.push_back(sId);
            }

            // If any, but not all, are on interior, project onto boundary
            if (phiX0 < 0.0 || phiX1 < 0.0 || phiX2 < 0.0) {
                if (phiX0 < 0) {
                    interpolateBoundaryLocation(x0, phiFun, nVals, bb);
                }
                if (phiX1 < 0) {
                    interpolateBoundaryLocation(x1, phiFun, nVals, bb);
                }
                if (phiX2 < 0) {
                    interpolateBoundaryLocation(x2, phiFun, nVals, bb);
                }
            }
        }

        // Any point epsilon close to the boundary is a boundary point
        const double eps = 1e-12;
        for (int i = 0; i < Vc->rows(); i++) {
            Eigen::Vector<double, D> x((*Vc)(i, Eigen::all));

            double phiVal = phiFun(x[0], x[1]);
            if (abs(phiVal) < eps) {
                boundaryMask->at(i) = bType;
            }
        }

    //     // Sort the removed Id's into reverse order for easy removal
    //     std::sort(idsToBeRemoved.begin(), idsToBeRemoved.end(), std::greater<int>());

    //     for (auto id = idsToBeRemoved.begin(); id != idsToBeRemoved.end(); ++id) {
    //         removeRow(*F, *id);
    //     }

    //     // Now, find all the used points.
    //     std::set<int> usedPnts{F->data(), F->data() + F->size()};

    //     Eigen::MatrixXd *VcTemp = new Eigen::MatrixXd(usedPnts.size(), D);
    //     int off = 0;

    //     // Build new set of points with offending elements removed.
    //     // Additionally, build map assigning new id's to the points in
    //     //  the simplex matrix.
    //     std::map<int, int> translator;
    //     for (int i = 0; i < Vc->size(); i++) {
    //         if (usedPnts.find(i) != usedPnts.end()) {
    //             (*VcTemp)(off, Eigen::all) = (*Vc)(i, Eigen::all);
    //             translator[i] = off;
    //             off++;
    //         }
    //     }

    //     // Make VC = VCtemp
    //     delete Vc;
    //     Vc = VcTemp;
    //     VcTemp = nullptr;

    //     // Iterate through simplex list, mapping the point ids
    //     for (int i = 0; i < F->rows(); i++) {
    //         for (int j = 0; j < D; j++) {
    //             (*F)(i, j) = translator[(*F)(i, j)];
    //         }
    //     }
    }
}

#endif