#include <iostream>
#include "./src/ADMMPG.h"
// #include "./src/Mesh2D.h"
#include "./src/PhaseM.h"
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include "./src/Mesh.h"
#include "./src/Assembly.h"
#include <cstdlib> 
#include <ctime> 
 
using namespace std;
#define D 2

void generateUniformRectMesh(unordered_map<string,double> params, Eigen::MatrixXd *Vc,
      Eigen::MatrixXi *F, Eigen::VectorXi *boundaryMask) {
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

    for (int i = 0; i < (nx+1)*(ny+1); i++) {
        int iOff = i % (nx+1);
        int jOff = i / (ny+1);
        bool boundaryPnt = (iOff == 0) || (iOff == nx) || (jOff == 0) || (jOff == ny);

        if (boundaryPnt) {
            (*boundaryMask)(i) = 1;
        }
    }
}

int main()
{
  srand(static_cast<unsigned int>(std::time(nullptr)));
  // Specify the monitor function
  // PhaseM<2> phaseM;
  PhaseM<2> *M = new PhaseM<2>();
  // MonitorFunction<2> *M = (MonitorFunction<2>*) new PhaseM<2>();
  // PhaseM M();

  // Parameters for the mesh
  std::unordered_map<std::string, double> params;
  int nx = 10;
  int ny = 10;
  int nPnts = (nx+1)*(ny+1) + nx*ny;
  params["nx"] = nx;
  params["ny"] = ny;
  params["nPnts"] = (params["nx"]+1)*(params["ny"]+1);
  params["d"] = D;
  double rho = 10;
  params["rho"] = rho;

  params["xa"] = 0.0;
  params["xb"] = 1.0;
  params["ya"] = 0.0;
  params["yb"] = 1.0;
  params["theta"] = 0.5;
  params["p"] = 1;
  double tau = 1e-1;
  params["tau"] = tau;

  Eigen::MatrixXd *Vc = nullptr;
  Eigen::MatrixXi *F = nullptr;
  Eigen::VectorXi *boundaryMask = nullptr;

  // Generate the initial mesh
  Vc = new Eigen::MatrixXd(nPnts, D);
  F = new Eigen::MatrixXi(4*nx*ny, D+1);

  cout << "Size of the matrix F " << F->rows() << endl;
//   assert(false);
  boundaryMask = new Eigen::VectorXi(Vc->rows());

  generateUniformRectMesh(params, Vc, F, boundaryMask);

  cout << "Creating the mesh" << endl;
  Mesh<D> adaptiveMesh(*Vc, *F, *boundaryMask, M, rho, tau);
  cout << "finished creating the mesh" << endl;

  // Create the solver
  double dt = 0.1;
  cout << "Creating the solver" << endl;
  ADMMPG<D> solver(dt, adaptiveMesh);
  cout << "FINISHED Creating the solver" << endl;

  clock_t start = clock();
  int nSteps = 5;
  cout << "Starting the time stepper" << endl;
  for (int i = 0; i < nSteps; i++) {
    solver.step(100, 1e-6);
    cout << "STEP = " << i << endl;
  }

  cout << "Took " << ((double)clock() - (double)start)
        / ((double)CLOCKS_PER_SEC) / ((double)nSteps) << " per steps" << endl;

  adaptiveMesh.outputPoints("points.txt");
  adaptiveMesh.outputSimplices("triangles.txt");

  delete M;
  delete Vc;
  delete F;
  delete boundaryMask;
}
