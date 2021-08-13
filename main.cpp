#include <iostream>
#include "./src/MeshIntegrator.h"
// #include "./src/Mesh2D.h"
#include "./src/PhaseM.h"
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include "./src/Mesh.h"
#include <cstdlib> 
#include <ctime> 
// #include <omp.h>
#include <unistd.h>
 
using namespace std;
#define D 2

void generateUniformRectMesh(unordered_map<string,double> params, Eigen::MatrixXd *Vc,
      Eigen::MatrixXi *F, vector<Mesh<2>::NodeType> *boundaryMask) {
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
            // boundaryMask->at(i) = Mesh<2>::NodeType::BOUNDARY_FIXED;
            boundaryMask->at(i) = Mesh<2>::NodeType::BOUNDARY_FREE;
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

int main() {
  srand(static_cast<unsigned int>(std::time(nullptr)));
  // Specify the monitor function
  PhaseM<2> *M = new PhaseM<2>();

  // Parameters for the mesh
  std::unordered_map<std::string, double> params;
  int nx = 20;
  int ny = 20;
  int nPnts = (nx+1)*(ny+1) + nx*ny;
  params["nx"] = nx;
  params["ny"] = ny;
  params["d"] = D;
  double rho = 5.0;
  params["rho"] = rho;

  params["xa"] = 0.0;
  params["xb"] = 1.0;
  params["ya"] = 0.0;
  params["yb"] = 1.0;
  params["theta"] = 0.5;
  params["p"] = 1;
  double tau = 1.0;
  params["tau"] = tau;

  Eigen::MatrixXd *Vc = nullptr;
  Eigen::MatrixXi *F = nullptr;
  vector<Mesh<2>::NodeType> *boundaryMask = nullptr;

  // Generate the initial mesh
  Vc = new Eigen::MatrixXd(nPnts, D);
  F = new Eigen::MatrixXi(4*nx*ny, D+1);

  boundaryMask = new vector<Mesh<2>::NodeType>(Vc->rows());

  generateUniformRectMesh(params, Vc, F, boundaryMask);

  cout << "Creating the mesh object" << endl;
  Mesh<D> adaptiveMesh(*Vc, *F, *boundaryMask, M, rho, tau);
  cout << "Finished Creating the mesh object" << endl;

  // Create the solver
  double dt = 0.1;
  cout << "Creating the sovler" << endl;
  MeshIntegrator<D> solver(dt, adaptiveMesh);
  cout << "FINISHED Creating the sovler" << endl;

  clock_t start = clock();
  int nSteps = 50; 
  double Ih;
  double Ihprev = INFINITY;
  int i;
  cout << "Running the solver" << endl;
  for (i = 0; i < nSteps; i++) {
    Ih = solver.step(100, 1e-5);
    cout << "Ih = " << Ih << endl;

    if (Ih >= Ihprev) {
        cout << "converged" << endl;
        break;
    }
    Ihprev = Ih;
  }

  nSteps = 1;
  cout << "Took " << ((double)clock() - (double)start)
        / ((double)CLOCKS_PER_SEC) << " seconds" << endl;
    cout << "Took " << i << " iters" << endl;

  adaptiveMesh.outputPoints("points.txt");
  adaptiveMesh.outputSimplices("triangles.txt");

  delete M;
  delete Vc;
  delete F;
  delete boundaryMask;
}
