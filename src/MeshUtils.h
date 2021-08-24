#ifndef MESH_UTILS_H
#define MESH_UTILS_H

#include <Eigen/Dense>
#include <vector>
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
        // x = (x - xMesh.at(0))/(xMesh.back() - xMesh.at(0));
        // y = (y - yMesh.at(0))/(yMesh.back() - yMesh.at(0));
        // cout << "norm = " << norm <<  endl;
        // cout << "x.size = " << xMesh.size() <<  endl;
        // cout << "(x, y) = (" << x << ", " << y << ")" << endl;
        coefs[0] = norm * (xMesh[1] - x)*(yMesh[1] - y);
        coefs[1] = norm * (x - xMesh[0])*(yMesh[1] - y);
        coefs[2] = norm * (xMesh[1] - x)*(y - yMesh[0]);
        coefs[3] = norm * (x - xMesh[0])*(y - yMesh[0]);
            
        // return (1/((xMesh[1]-xMesh[0])*(yMesh[1] - yMesh[0])))*(
        //             f[0]*(xMesh[1] - x)*(yMesh[1] - y)
        //                 + f[1]*(x - xMesh[0])*(yMesh[1] - y)
        //                 + f[2]*(xMesh[1] - x)*(y - yMesh[0])
        //                 + f[3]*(x - xMesh[0])*(y - yMesh[0]) );
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
        // cout << "numBound = " << numBound << endl;

        // // Now that we have the points on a grid, map them to a circle
        // for (int i = 0; i < Vc->rows(); i++) {
        //     Eigen::Vector<double, D> temp;// = (*Vc)(i, Eigen::all);

        //     temp(0) = -1 + 2*(*Vc)(i, 0);
        //     temp(1) = -1 + 2*(*Vc)(i, 1);

        //     (*Vc)(i, Eigen::all) = temp;
        // }

        // for (int i = 0; i < Vc->rows(); i++) {
        //     Eigen::Vector<double, D> temp;// = (*Vc)(i, Eigen::all);

        //     temp(0) = ((*Vc)(i, 0))*sqrt(1 - pow(((*Vc)(i, 1)), 2.0)/2.0);
        //     temp(1) = ((*Vc)(i, 1))*sqrt(1 - pow(((*Vc)(i, 0)), 2.0)/2.0);

        //     (*Vc)(i, Eigen::all) = temp;
        // }

        // for (int i = 0; i < Vc->rows(); i++) {
        //     Eigen::Vector<double, D> temp;// = (*Vc)(i, Eigen::all);

        //     temp(0) = (1 + (*Vc)(i, 0))/2.0;
        //     temp(1) = (1 + (*Vc)(i, 1))/2.0;

        //     (*Vc)(i, Eigen::all) = temp;
        // }

    }
}

#endif