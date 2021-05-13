#include <iostream>
#include "./src/MeshIntegrator.h"
#include "./src/PhaseM.h"
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include "./src/Mesh.h"
#include <cstdlib> 
#include <ctime> 
#include <unistd.h>
 
using namespace std;
#define D 3

void generateUniformRectMesh(unordered_map<string,double> params, Eigen::MatrixXd *Vc,
      Eigen::MatrixXi *F, vector<Mesh<3>::NodeType> *boundaryMask) {
    int nx = (int) params["nx"];
    int ny = (int) params["ny"];
    int nz = (int) params["nz"];

    int xa = params["xa"];
    int xb = params["xb"];
    int ya = params["ya"];
    int yb = params["yb"];
    int za = params["za"];
    int zb = params["zb"];

    double hx = (xb - xa)/((double)nx);
    double hy = (yb - ya)/((double)ny);
    double hz = (zb - za)/((double)nz);

    // Append the cube vertices
    int off = 0;
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
        boundaryMask->at(i) = Mesh<3>::NodeType::INTERIOR;
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
                boundaryMask->at(off) = Mesh<3>::NodeType::BOUNDARY_FREE;
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
                boundaryMask->at(off) = Mesh<3>::NodeType::BOUNDARY_FIXED;
            }
        }
    }
}

int main()
{
  srand(static_cast<unsigned int>(std::time(nullptr)));
  // Specify the monitor function
  PhaseM<D> *M = new PhaseM<D>();

  // Parameters for the mesh
  std::unordered_map<std::string, double> params;
  int nx = 20;
  int ny = 20;
  int nz = 20;
  int nPnts = (nx+1)*(ny+1)*(nz+1) + nx*ny*nz;
  params["nx"] = nx;
  params["ny"] = ny;
  params["nz"] = nz;
  params["nPnts"] = (params["nx"]+1)*(params["ny"]+1)*(params["nz"]+1);
  params["d"] = D;
  double rho = 10.0;
  params["rho"] = rho;

  params["xa"] = 0.0;
  params["xb"] = 1.0;
  params["ya"] = 0.0;
  params["yb"] = 1.0;
  params["za"] = 0.0;
  params["zb"] = 1.0;
  params["theta"] = 0.5;
  params["p"] = 1;
  double tau = 1e-2;
  params["tau"] = tau;

  Eigen::MatrixXd *Vc = nullptr;
  Eigen::MatrixXi *F = nullptr;
  vector<Mesh<3>::NodeType> *boundaryMask = nullptr;

  // Generate the initial mesh
  Vc = new Eigen::MatrixXd(nPnts, D);
  F = new Eigen::MatrixXi(12*nx*ny*nz, D+1);

  cout << "Size of the matrix F " << F->rows() << endl;
//   assert(false);
  boundaryMask = new vector<Mesh<3>::NodeType>(Vc->rows());

  cout << "gen init mesh" << F->rows() << endl;
  generateUniformRectMesh(params, Vc, F, boundaryMask);
  cout << "FINISHED gen init mesh" << F->rows() << endl;
  
  cout << "Creating the mesh" << endl;
  Mesh<D> adaptiveMesh(*Vc, *F, *boundaryMask, M, rho, tau);
  cout << "finished creating the mesh" << endl;

  // Create the solver
  double dt = 0.05;
  cout << "Creating the solver" << endl;
  MeshIntegrator<D> solver(dt, adaptiveMesh);
  cout << "FINISHED Creating the solver" << endl;

  clock_t start = clock();
  int nSteps = 50; 
  cout << "Starting the time stepper" << endl;
  double Ih;
  double Ihprev = INFINITY;
  for (int i = 0; i < nSteps; i++) {
    Ih = solver.step(100, 1e-3);
    cout << "STEP = " << i << endl;
    cout << "Ih = " << Ih << endl;

    if (Ih >= Ihprev) {
        cout << "converged" << endl;
        break;
    }
    Ihprev = Ih;
  }

//   nSteps = 1;
//   cout << "Took " << ((double)clock() - (double)start)
//         / ((double)CLOCKS_PER_SEC) / ((double)nSteps) << " per steps" << endl;

  adaptiveMesh.outputPoints("points3d.txt");
  adaptiveMesh.outputSimplices("tets.txt");
  adaptiveMesh.outputBoundaryNodes("boundaryPnts.txt");

  delete M;
  delete Vc;
  delete F;
  delete boundaryMask;
}
