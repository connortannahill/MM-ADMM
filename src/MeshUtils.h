#ifndef MESH_UTILS_H
#define MESH_UTILS_H

#include <Eigen/Dense>
#include <vector>

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
        int guess =  (int)((w - w_mesh.at(0))/(w_mesh.at(1) - w_mesh.at(0)));
        if (guess < 0) {
            guess = 0;
        } else if (guess > (int)((w_mesh.back() - w_mesh.at(0))/(w_mesh.at(1) - w_mesh.at(0)))) {
            guess = (int)((w_mesh.back() - w_mesh.at(0))/(w_mesh.at(1) - w_mesh.at(0)));
        }
        return guess;
    }

    inline void biLinearInterpolation(double x, double y, vector<double> &xMesh,
                                            vector<double> &yMesh, double *coefs) {
        coefs[0] = (1/((xMesh.at(1)-xMesh.at(0))*(yMesh.at(1) - yMesh.at(0)))) * (xMesh.at(1) - x)*(yMesh.at(1) - y);
        coefs[1] = (1/((xMesh.at(1)-xMesh.at(0))*(yMesh.at(1) - yMesh.at(0)))) * (x - xMesh.at(0))*(yMesh.at(1) - y);
        coefs[2] = (1/((xMesh.at(1)-xMesh.at(0))*(yMesh.at(1) - yMesh.at(0)))) * (xMesh.at(1) - x)*(y - yMesh.at(0));
        coefs[3] = (1/((xMesh.at(1)-xMesh.at(0))*(yMesh.at(1) - yMesh.at(0)))) * (x - xMesh.at(0))*(y - yMesh.at(0));
            
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
}

#endif