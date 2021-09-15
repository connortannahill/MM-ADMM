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

double circlePhi(double x, double y) {
  double r = 0.4;
  double cx = 0.5;
  double cy = 0.5;
  double xval = (x - cx);
  double yval = (y - cy);
  return sqrt(xval*xval + yval*yval) - r;
}

void setUpLevelSetExperiment(string testName, ifstream &inputFile,
    int numThreads, MonitorFunction<D> *mon) {
  int nx, ny, nSteps;
  double dt, tau, rho;
  double xa, xb, ya, yb;

  int boundaryType;
  inputFile >> boundaryType;

  NodeType type;
  if (boundaryType == 0) {
    type = NodeType::BOUNDARY_FREE;
  } else {
    type = NodeType::BOUNDARY_FIXED;
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

  cout << "done params" << endl;

  Eigen::MatrixXd *Vc = new Eigen::MatrixXd(1, D);
  Eigen::MatrixXi *F = new Eigen::MatrixXi(1, D+1);
  vector<NodeType> *boundaryMask = new vector<NodeType>();

  // Put all
  std::function<double(double,double)> phi = circlePhi;
  std::vector<int> nVals;
  nVals.push_back(nx);
  nVals.push_back(ny);

  std::vector<std::tuple<double, double>> bb;
  bb.push_back(std::make_tuple(xa, xb));
  bb.push_back(std::make_tuple(ya, yb));

  cout << "done making the weird stuff" << endl;

  utils::meshFromLevelSetFun(phi, nVals, bb, Vc, F, boundaryMask, type);

  cout << "done making mesh from the level set fcn" << endl;

  cout << endl;
  cout << numThreads << endl;

  Mesh<D> adaptiveMesh(*Vc, *F, *boundaryMask, mon,
      numThreads, rho, tau);
  
  cout << "FInished making the mesh" << endl;

  /** Ouptut the mesh */
  string outDir = "./Experiments/Results/" + testName;
  int stat = system(("mkdir -p " + outDir).c_str());
  if (stat != 0) {
    cout << "mkdir failed" << endl;
    assert(false);
  }

  string pointsOutDir = outDir + "/points.txt";
  string triangleOutDir = outDir + "/triangles.txt";
  adaptiveMesh.outputPoints(pointsOutDir.c_str());
  adaptiveMesh.outputSimplices(triangleOutDir.c_str());

  std::ofstream outFile;
  string phiOutDir = outDir + "/phi.txt";
  outFile.open(phiOutDir.c_str());

  std::vector<double> x;
  std::vector<double> y;

  cout << "linspace" << endl;
  utils::linspace(xa, xb, nx, x);
  utils::linspace(ya, yb, ny, y);
  cout << "DONE linspace" << endl;

  cout << "outting" << endl;
  for (int i = 0; i < x.size(); i++) {
    for (int j = 0; j < y.size(); j++) {
      outFile << x.at(i) << ", " << y.at(j) << ", " << phi(x.at(i), y.at(j)) << endl;
    }
  }
  cout << "fINISHED outting" << endl;
  outFile.close();

  delete Vc;
  delete F;
  delete boundaryMask;

  cout << "Done delete" << endl;
}

void setUpBoxExperiment(string testName, ifstream &inputFile,
    int numThreads, MonitorFunction<D> *mon) {
  int nx, ny, nSteps;
  double dt, tau, rho;
  double xa, xb, ya, yb;

  int boundaryType;
  inputFile >> boundaryType;

  NodeType type;
  if (boundaryType == 0) {
    type = NodeType::BOUNDARY_FREE;
  } else {
    type = NodeType::BOUNDARY_FIXED;
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
  vector<NodeType> *boundaryMask = nullptr;

  // Generate the initial mesh
  Vc = new Eigen::MatrixXd(nPnts, D);
  F = new Eigen::MatrixXi(4*nx*ny, D+1);

  boundaryMask = new vector<NodeType>(Vc->rows());

  utils::generateUniformRectMesh(params, Vc, F, boundaryMask,
      type);

  Mesh<D> adaptiveMesh(*Vc, *F, *boundaryMask, mon,
      numThreads, rho, tau);

  // Create the solver
  MeshIntegrator<D> solver(dt, adaptiveMesh);

  clock_t start = clock();
  double Ih;
  double Ihprev = INFINITY;
  int i;
  for (i = 0; i < nSteps; i++) {
    Ih = solver.step(50, 1e-4);

    // cout << "Ih = " << Ih << endl;
    // cout << "IhDiff = " << Ih - Ihprev << endl;
    double dIdt = (Ih - Ihprev) / dt;
    cout << "d/dt = " << (Ih - Ihprev) / dt << endl;
    cout << "Ih = " << Ih << endl;

    // if (i != 0 && Ih >= Ihprev) {
    if (i != 0 && (abs(dIdt) < 1e-3 || dIdt > 0)) {
        cout << "converged" << endl;
        break;
    }
    Ihprev = Ih;
  }

  cout << "Took " << ((double)clock() - (double)start)
        / ((double)CLOCKS_PER_SEC) << " seconds" << endl;
  cout << "Took " << i << " iters" << endl;

  string outDir = "./Experiments/Results/" + testName;
  int stat = system(("mkdir -p " + outDir).c_str());
  if (stat != 0) {
    cout << "mkdir failed" << endl;
    assert(false);
  }

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
  assert(argc <= 3);
  string inFileName = argv[1];
  cout << inFileName << endl;

  int numThreads = 1;
  if (argc == 3) {
    numThreads = atoi(argv[2]);
  }

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
    setUpBoxExperiment(inFileName, inputFile, numThreads, Mvals.at(monType-1));
  } else if (testType.compare("LevelSet") == 0) {
    setUpLevelSetExperiment(inFileName, inputFile, numThreads, Mvals.at(monType-1));
  } else {
    assert(false);
  }

  inputFile.close();
}
