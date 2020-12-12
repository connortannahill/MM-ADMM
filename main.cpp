#include <iostream>
#include "./src/ADMMPG.h"
#include "./src/Mesh2D.h"
#include "./src/PhaseM.h"
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
 
using namespace std;
 
int main()
{
  // Specify the monitor function
  PhaseM *M = new PhaseM();
  // PhaseM M();

  // Parameters for the mesh
  std::unordered_map<std::string, double> params;
  params["nx"] = 10;
  params["ny"] = 10;
  params["nPnts"] = (params["nx"]+1)*(params["ny"]+1);
  params["d"] = 2;
  params["rho"] = 10;
  // params["rho"] = 0.0;

  params["xa"] = 0.0;
  params["xb"] = 1.0;
  params["ya"] = 0.0;
  params["yb"] = 1.0;
  params["theta"] = 0.5;
  params["p"] = 1;
  params["tau"] = 1e-2;
  // params["tau"] = 1e-1;

  // Make a square mesh
  Mesh2D adaptiveMesh(M, params);

  // for (int i = 0; i < 300; i++) {
  //     adaptiveMesh.eulerStep(0.001);
  // }

  // Create the solver
  double dt = 0.01;
  ADMMPG solver(dt, adaptiveMesh);


  clock_t start = clock();
  int nSteps = 30;
  for (int i = 0; i < nSteps; i++) {
    solver.step(200, 1e-5);
  }

  cout << "Time per step = " << ((clock() - start)/((double)CLOCKS_PER_SEC))/((double)nSteps);

  delete M;

  adaptiveMesh.outputPoints("points.txt");
  adaptiveMesh.outputTriangles("triangles.txt");
}
