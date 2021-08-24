#include <iostream>
#include "./src/MeshIntegrator.h"
#include "./src/Mesh.h"
#include <unordered_map>
#include "Experiments/TestMonitors/MEx1.h"
#include "Experiments/TestMonitors/MEx2.h"
#include "Experiments/TestMonitors/MEx3.h"
#include <string>
#include "./src/Mesh.h"
#include "./src/MonitorFunction.h"
#include <cstdlib> 
#include <ctime> 
#include <fstream>
#include <unistd.h>
#include "./src/MeshUtils.h"
 
using namespace std;
#define D 2

void setUpBoxExperiment(string testName, ifstream &inputFile,
    MonitorFunction<D> *mon) {
  int nx, ny, nSteps;
  double dt, tau, rho;
  double xa, xb, ya, yb;

  int boundaryType
  inputFile >> boundaryType;

  Mesh<D>::NodeType type;
  if (boundaryType == 0) {
    type = Mesh<D>::NodeType::BOUNDARY_FREE;
  } else {
    type = Mesh<D>::NodeType::BOUNDARY_FIXED;
  }

  inputFile >> nSteps;
  inputFile >> dt;
  inputFile >> tau;
  inputFile >> rho;

  int nPnts;
  inputFile >> nx >> ny;
  nPnts = (nx+1)*(ny+1) + nx*ny;

  inputFile >> xa >> xb >> ya >> yb;

  // Parameters for the mesh
  std::unordered_map<std::string, double> params;
  params["nx"] = nx;
  params["ny"] = ny;
  params["d"] = D;
  params["rho"] = rho;
  params["tau"] = tau;

  params["xa"] = xa;
  params["xb"] = xb;
  params["ya"] = ya;
  params["yb"] = yb;

  Eigen::MatrixXd *Vc = nullptr;
  Eigen::MatrixXi *F = nullptr;
  vector<Mesh<2>::NodeType> *boundaryMask = nullptr;

  // Generate the initial mesh
  Vc = new Eigen::MatrixXd(nPnts, D);
  F = new Eigen::MatrixXi(4*nx*ny, D+1);

  boundaryMask = new vector<Mesh<2>::NodeType>(Vc->rows());

  utils::generateUniformRectMesh(params, Vc, F, boundaryMask,
      type);

  Mesh<D> adaptiveMesh(*Vc, *F, *boundaryMask, mon, rho, tau);

  // Create the solver
  MeshIntegrator<D> solver(dt, adaptiveMesh);

  clock_t start = clock();
  double Ih;
  double Ihprev = INFINITY;
  int i;
  for (i = 0; i < nSteps; i++) {
    Ih = solver.step(100, 1e-4);


    if (Ih >= Ihprev) {
        cout << "converged" << endl;
        break;
    }
    cout << "Ih = " << Ih << endl;
    cout << "IhDiff = " << abs(Ih - Ihprev) << endl;
    Ihprev = Ih;
  }

  cout << "Took " << ((double)clock() - (double)start)
        / ((double)CLOCKS_PER_SEC) << " seconds" << endl;
  cout << "Took " << i << " iters" << endl;

  string outDir = "./Experiments/Results/" + testName;
  system(("mkdir -p " + outDir).c_str());

  string pointsOutDir = outDir + "/points.txt";
  string triangleOutDir = outDir + "/triangles.txt";
  adaptiveMesh.outputPoints(pointsOutDir.c_str());
  adaptiveMesh.outputSimplices(triangleOutDir.c_str());

  delete Vc;
  delete F;
  delete boundaryMask;
}

int main(int argc, char *argv[]) {
  srand(static_cast<unsigned int>(std::time(nullptr)));

  // Read in the input file
  assert(argc == 2);
  string inFileName = argv[1];
  cout << inFileName << endl;

  ifstream inputFile("Experiments/InputFiles/" + inFileName);

  string testType;
  getline(inputFile, testType);

  // Specify the monitor function
  auto *M1 = new MEx1<D>();
  auto *M2 = new MEx2<D>();
  auto *M3 = new MEx3<D>();
  vector<MonitorFunction<D>*> Mvals;
  Mvals.push_back(M1);
  Mvals.push_back(M2);
  Mvals.push_back(M3);

  // Get monitor id
  int monType;
  inputFile >> monType;
  cout << "mon type = " << monType << endl;

  if (testType.compare("SquareGrid") == 0) {
    setUpBoxExperiment(inFileName, inputFile, Mvals.at(monType-1));
  } else {
    assert(false);
  }

  inputFile.close();
  




//   // Parameters for the mesh
  // std::unordered_map<std::string, double> params;
//   int nx = 20;
//   int ny = 20;
//   int nPnts = (nx+1)*(ny+1) + nx*ny;
//   params["nx"] = nx;
//   params["ny"] = ny;
//   params["d"] = D;
//   double rho = 10;
//   params["rho"] = rho;

//   params["xa"] = 0.0;
//   params["xb"] = 1.0;
//   params["ya"] = 0.0;
//   params["yb"] = 1.0;
//   params["tau"] = tau;

//   Eigen::MatrixXd *Vc = nullptr;
//   Eigen::MatrixXi *F = nullptr;
//   vector<Mesh<2>::NodeType> *boundaryMask = nullptr;

//   // Generate the initial mesh
//   Vc = new Eigen::MatrixXd(nPnts, D);
//   F = new Eigen::MatrixXi(4*nx*ny, D+1);

//   boundaryMask = new vector<Mesh<2>::NodeType>(Vc->rows());

//   utils::generateUniformRectMesh(params, Vc, F, boundaryMask, Mesh<2>::NodeType::BOUNDARY_FREE);

//   cout << "Creating the mesh object" << endl;
//   Mesh<D> adaptiveMesh(*Vc, *F, *boundaryMask, M, rho, tau);
//   cout << "Finished Creating the mesh object" << endl;

//   // Create the solver
//   double dt = 0.1;
//   cout << "Creating the sovler" << endl;
//   MeshIntegrator<D> solver(dt, adaptiveMesh);
//   cout << "FINISHED Creating the sovler" << endl;

//   clock_t start = clock();
//   int nSteps = 50; 
//   double Ih;
//   double Ihprev = INFINITY;
//   int i;
//   cout << "Running the solver" << endl;
//   for (i = 0; i < nSteps; i++) {
//     Ih = solver.step(100, 1e-6);
//     cout << "Ih = " << Ih << endl;

//     if (Ih >= Ihprev) {
//         cout << "converged" << endl;
//         break;
//     }
//     Ihprev = Ih;
//   }

//   nSteps = 1;
//   cout << "Took " << ((double)clock() - (double)start)
//         / ((double)CLOCKS_PER_SEC) << " seconds" << endl;
//     cout << "Took " << i << " iters" << endl;

//   adaptiveMesh.outputPoints("points.txt");
//   adaptiveMesh.outputSimplices("triangles.txt");

//   delete M;
//   delete Vc;
//   delete F;
//   delete boundaryMask;
}
