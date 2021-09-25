#include <iostream>
#include "./src/MeshIntegrator.h"
#include "./src/Mesh.h"
#include <unordered_map>
#include "Experiments/TestMonitors/MEx0.h"
#include "Experiments/TestMonitors/MEx1.h"
#include "Experiments/TestMonitors/MEx2.h"
#include "Experiments/TestMonitors/MEx3.h"
#include "Experiments/TestMonitors/MEx13D.h"
#include "Experiments/TestMonitors/MEx23D.h"
#include "Experiments/TestMonitors/MEx33D.h"
#include <string>
#include "./src/Mesh.h"
#include "./src/MonitorFunction.h"
#include <cstdlib> 
#include <ctime> 
#include <fstream>
#include <unistd.h>
#include "./src/MeshUtils.h"

using namespace std;

double circlePhi(double x, double y) {
  double r = 0.4;
  double cx = 0.5;
  double cy = 0.5;
  double xval = (x - cx);
  double yval = (y - cy);
  // return sqrt(xval*xval + yval*yval) - r;
  return xval*xval + yval*yval - r*r;
}

double spherePhi(double x, double y, double z) {
  double r = 0.4;
  double cx = 0.5;
  double cy = 0.5;
  double cz = 0.5;
  double xval = (x - cx);
  double yval = (y - cy);
  double zval = (z - cz);
  // return sqrt(xval*xval + yval*yval) - r;
  return xval*xval + yval*yval + zval*zval - r*r;
}

double heartPhi(double x, double y) {
  // double cx = 0.5;
  // double cy = 0.4;
  double cx = 0.5;
  double cy = 2.4;
  double xval = (x - cx);
  double yval = (y - cy);

  // return pow(pow(xval, 2.0) + pow(yval, 2.0) - 0.1, 3.0) - pow(xval, 2.0)*pow(yval, 3.0);
  return pow(yval - (2.0*(abs(xval)+pow(xval, 2)-6))/(3.0*(abs(xval)+pow(xval, 2)+2)), 2.0) + pow(xval, 2) - 0.1;
}

// Doesnt work
double shoulderPhi(double x, double y) {
  double n = 500;
  double r1 = 0.4;
  double cx1 = 0.5;
  double cy1 = 0.5;
  double xval1 = (x - cx1);
  double yval1 = (y - cy1);

  double phi1 = pow(xval1, n) + pow(yval1, n) - pow(r1, n);

  double r2 = 0.2;
  double cx2 = 0.675;
  double cy2 = 0.675;
  double xval2 = (x - cx2);
  double yval2 = (y - cy2);

  double phi2 = pow(xval2, n) + pow(yval2, n) - pow(r2, n);

  return max(phi1, phi2);
  // return phi2;
}

template <typename T>
void outputVecToFile(string fileName, vector<T> &outVec) {
  std::ofstream outFile;
  outFile.open(fileName);

  for (auto val = outVec.begin(); val != outVec.end(); ++val) {
    outFile << *val << endl;
  }
  outFile.close();
}

template <int D>
void setUpLevelSetExperiment(string testName, ifstream &inputFile,
    int numThreads, MonitorFunction<D> *mon) {
  int nx, ny, nz, nSteps;
  double dt, tau, rho;
  double xa, xb, ya, yb, za, zb;

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

  if (D == 2) {
    inputFile >> nx >> ny;
    inputFile >> xa >> xb >> ya >> yb;
    nz = 0;
    za = 0;
    zb = 0;
  } else {
    inputFile >> nx >> ny >> nz;
    inputFile >> xa >> xb >> ya >> yb >> za >> zb;
  }

  // Parameters for the mesh
  std::unordered_map<std::string, double> params;
  params["nx"] = nx;
  params["ny"] = ny;
  params["nz"] = nz;
  params["d"] = D;
  params["rho"] = rho;
  params["tau"] = tau;

  params["xa"] = xa;
  params["xb"] = xb;
  params["ya"] = ya;
  params["yb"] = yb;
  params["za"] = za;
  params["zb"] = zb;

  Eigen::MatrixXd *Vc = new Eigen::MatrixXd(1, D);
  Eigen::MatrixXd *Vp = new Eigen::MatrixXd(1, D);
  Eigen::MatrixXi *F = new Eigen::MatrixXi(1, D+1);
  vector<NodeType> *boundaryMask = new vector<NodeType>();

  string outDir = "./Experiments/Results/" + testName;

  if (D == 2) {
    std::function<double(double,double)> phi = circlePhi;

    std::vector<int> nVals;
    nVals.push_back(nx);
    nVals.push_back(ny);

    std::vector<std::tuple<double, double>> bb;
    bb.push_back(std::make_tuple(xa, xb));
    bb.push_back(std::make_tuple(ya, yb));

    utils::meshFromLevelSetFun(phi, nVals, bb, Vc, Vp, F, boundaryMask, type);

    std::ofstream outFile;
    string phiOutDir = outDir + "/phi.txt";
    outFile.open(phiOutDir.c_str());

    std::vector<double> x;
    std::vector<double> y;

    utils::linspace(xa, xb, nx, x);
    utils::linspace(ya, yb, ny, y);

    for (uint32_t i = 0; i < x.size(); i++) {
      for (uint32_t j = 0; j < y.size(); j++) {
        outFile << x.at(i) << ", " << y.at(j) << ", " << phi(x.at(i), y.at(j)) << endl;
      }
    }
    outFile.close();
  } else {
    std::function<double(double,double,double)> phi = spherePhi;

    std::vector<int> nVals;
    nVals.push_back(nx);
    nVals.push_back(ny);
    nVals.push_back(nz);

    std::vector<std::tuple<double, double>> bb;
    bb.push_back(std::make_tuple(xa, xb));
    bb.push_back(std::make_tuple(ya, yb));
    bb.push_back(std::make_tuple(za, zb));

    utils::meshFromLevelSetFun(phi, nVals, bb, Vc, Vp, F, boundaryMask, type);

    std::ofstream outFile;
    string phiOutDir = outDir + "/phi.txt";
    outFile.open(phiOutDir.c_str());

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    utils::linspace(xa, xb, nx, x);
    utils::linspace(ya, yb, ny, y);
    utils::linspace(za, zb, nz, z);

    for (uint32_t k = 0; k < z.size(); k++) {
      for (uint32_t j = 0; j < y.size(); j++) {
        for (uint32_t i = 0; i < x.size(); i++) {
          outFile << x.at(i) << ", " << y.at(j) << ", " << z.at(k) << ", " << phi(x.at(i), y.at(j), z.at(k)) << endl;
        }
      }
    }
    outFile.close();
  }

  Mesh<D> adaptiveMesh(*Vp, *F, *boundaryMask, mon,
      numThreads, rho, tau);

  // Create and run the solver
  MeshIntegrator<D> solver(dt, adaptiveMesh);

  std::vector<double> Ivals;

  clock_t start = clock();
  double Ih;
  double Ihprev = INFINITY;
  int i;
  for (i = 0; i < nSteps; i++) {
    Ih = solver.step(50, 1e-4);

    Ivals.push_back(Ih);

    double dIdt = (Ih - Ihprev) / dt;
    cout << "d/dt = " << (Ih - Ihprev) / dt << endl;
    cout << "Ih = " << Ih << endl;

    if (i != 0 && (abs(dIdt) < 1e-4 || dIdt > 0)) {
        cout << "converged" << endl;
        break;
    }
    Ihprev = Ih;
  }

  cout << "Took " << ((double)clock() - (double)start)
        / ((double)CLOCKS_PER_SEC) << " seconds" << endl;
  cout << "Took " << i << " iters" << endl;
  
  /** Output the mesh */
  int stat = system(("mkdir -p " + outDir).c_str());
  if (stat != 0) {
    cout << "mkdir failed" << endl;
    assert(false);
  }

  string pointsOutDir = outDir + "/points.txt";
  string triangleOutDir = outDir + "/triangles.txt";
  string IhDir = outDir + "Ih.txt";
  adaptiveMesh.outputPoints(pointsOutDir.c_str());
  adaptiveMesh.outputSimplices(triangleOutDir.c_str());
  outputVecToFile(IhDir, Ivals);

  delete Vc;
  delete F;
  delete boundaryMask;
}

template <int D>
void setUpShoulderExperiment(string testName, ifstream &inputFile,
    int numThreads, MonitorFunction<D> *mon) {
  int nx, ny, nz, nSteps;
  double dt, tau, rho;
  double xa, xb, ya, yb, za, zb;

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
  if (D == 2) {
    inputFile >> nx >> ny;
    nPnts = (nx+1)*(ny+1) + nx*ny;
    inputFile >> xa >> xb >> ya >> yb;

    za = 0;
    zb = 0;
    nz = 0;
  } else {
    inputFile >> nx >> ny >> nz;
    nPnts = (nx+1)*(ny+1)*(nz+1) + nx*ny*nz;
    inputFile >> xa >> xb >> ya >> yb >> za >> zb;
  }


  // Parameters for the mesh
  std::unordered_map<std::string, double> params;
  params["nx"] = nx;
  params["ny"] = ny;
  params["nz"] = nz;
  params["d"] = D;
  params["rho"] = rho;
  params["tau"] = tau;

  params["xa"] = xa;
  params["xb"] = xb;
  params["ya"] = ya;
  params["yb"] = yb;
  params["za"] = za;
  params["zb"] = zb;

  Eigen::MatrixXd *Vc = nullptr;
  Eigen::MatrixXd *Vp = nullptr;
  Eigen::MatrixXi *F = nullptr;
  vector<NodeType> *boundaryMask = nullptr;

  // Generate the initial mesh
  Vp = new Eigen::MatrixXd(nPnts, D);
  Vc = new Eigen::MatrixXd(nPnts, D);
  if (D == 2) {
    F = new Eigen::MatrixXi(4*nx*ny, D+1);
  } else {
    F = new Eigen::MatrixXi(12*nx*ny*nz, D+1);
  }

  boundaryMask = new vector<NodeType>(Vc->rows());

  utils::generateUniformRectMesh<D>(params, Vc, F, boundaryMask,
      type);
  
  *Vp = *Vc;
    
  // Now we remove everything in the top-right corner
  double cx = (xa + xb)/2.0;
  double cy = (ya + yb)/2.0;
  double cz = (za + zb)/2.0;
  std::vector<int> idsToBeRemoved;

  // Any simplex with centroid in the shoulder region are removed
  const double EPS = 1e-16;
  for (int i = 0; i < F->rows(); i++) {
    Eigen::Vector<double, D> x0((*Vc)((*F)(i, 0), Eigen::all));
    Eigen::Vector<double, D> x1((*Vc)((*F)(i, 1), Eigen::all));
    Eigen::Vector<double, D> x2((*Vc)((*F)(i, 2), Eigen::all));
    Eigen::Vector<double, D> x3;
    if (D == 3) {
      x3 = (*Vc)((*F)(i, 3), Eigen::all);
    }
    
    Eigen::Vector<double, D> c;
    if (D == 2) {
      c = (1.0/3.0)*(x0 + x1 + x2);
    } else {
      c = (1.0/4.0)*(x0 + x1 + x2 + x3);
    }

    if (D == 2) {
      if (c[0] > cx && c[1] > cy) {
        idsToBeRemoved.push_back(i);

        // Label all of these as boundary points, as it wont matter.
        // In the case of equality, fix the corner points
        bool cond1, cond2, cond3;
        cond1 = abs(x0[0] - cx) < EPS && abs(x0[1] - cy) < EPS;
        cond2 = abs(x0[0] - cx) < EPS && abs(x0[1] - yb) < EPS;
        cond3 = abs(x0[0] - xb) < EPS && abs(x0[1] - cy) < EPS;
        boundaryMask->at((*F)(i, 0)) = (cond1 || cond2 || cond3) ? NodeType::BOUNDARY_FIXED : type;
        bool c1 = (cond1 || cond2 || cond3);

        cond1 = abs(x1[0] - cx) < EPS && abs(x1[1] - cy) < EPS;
        cond2 = abs(x1[0] - cx) < EPS && abs(x1[1] - yb) < EPS;
        cond3 = abs(x1[0] - xb) < EPS && abs(x1[1] - cy) < EPS;
        boundaryMask->at((*F)(i, 1)) = (cond1 || cond2 || cond3) ? NodeType::BOUNDARY_FIXED : type;
        bool c2 = (cond1 || cond2 || cond3);

        cond1 = abs(x2[0] - cx) < EPS && abs(x2[1] - cy) < EPS;
        cond2 = abs(x2[0] - cx) < EPS && abs(x2[1] - yb) < EPS;
        cond3 = abs(x2[0] - xb) < EPS && abs(x2[1] - cy) < EPS;
        boundaryMask->at((*F)(i, 2)) = (cond1 || cond2 || cond3) ? NodeType::BOUNDARY_FIXED : type;
        bool c3 = (cond1 || cond2 || cond3);

        assert(c1 + c2 + c3 <= 1);
      }
    } else {
      if (c[0] > cx && c[1] > cy && c[2] > cz) {
        idsToBeRemoved.push_back(i);

        // Label all of these as boundary points, as it wont matter.
        // In the case of equality, fix the corner points
        bool cond1, cond2, cond3, cond4, cond5, cond6, cond7;
        cond1 = abs(x0[0] - cx) < EPS && abs(x0[2] - cz) < EPS;
        cond2 = abs(x0[0] - cx) < EPS && abs(x0[2] - zb) < EPS;
        cond3 = abs(x0[0] - xb) < EPS && abs(x0[2] - cz) < EPS;
        cond4 = abs(x0[1] - ya) < EPS && abs(x0[2] - cz) < EPS;
        cond5 = abs(x0[1] - yb) < EPS && abs(x0[2] - cz) < EPS;
        cond6 = abs(x0[0] - cx) < EPS && abs(x0[1] - ya) < EPS;
        cond7 = abs(x0[0] - cx) < EPS && abs(x0[1] - yb) < EPS;
        bool c0 = (cond1 || cond2 || cond3 || cond4 || cond5 ||
                    cond6 || cond7);
        boundaryMask->at((*F)(i, 0)) = c0 ? NodeType::BOUNDARY_FIXED : type;

        cond1 = abs(x1[0] - cx) < EPS && abs(x1[2] - cz) < EPS;
        cond2 = abs(x1[0] - cx) < EPS && abs(x1[2] - zb) < EPS;
        cond3 = abs(x1[0] - xb) < EPS && abs(x1[2] - cz) < EPS;
        cond4 = abs(x1[1] - ya) < EPS && abs(x1[2] - cz) < EPS;
        cond5 = abs(x1[1] - yb) < EPS && abs(x1[2] - cz) < EPS;
        cond6 = abs(x1[0] - cx) < EPS && abs(x1[1] - ya) < EPS;
        cond7 = abs(x1[0] - cx) < EPS && abs(x1[1] - yb) < EPS;
        bool c1 = (cond1 || cond2 || cond3 || cond4 || cond5 ||
                    cond6 || cond7);
        boundaryMask->at((*F)(i, 1)) = c1 ? NodeType::BOUNDARY_FIXED : type;

        cond1 = abs(x2[0] - cx) < EPS && abs(x2[2] - cz) < EPS;
        cond2 = abs(x2[0] - cx) < EPS && abs(x2[2] - zb) < EPS;
        cond3 = abs(x2[0] - xb) < EPS && abs(x2[2] - cz) < EPS;
        cond4 = abs(x2[1] - ya) < EPS && abs(x2[2] - cz) < EPS;
        cond5 = abs(x2[1] - yb) < EPS && abs(x2[2] - cz) < EPS;
        cond6 = abs(x2[0] - cx) < EPS && abs(x2[1] - ya) < EPS;
        cond7 = abs(x2[0] - cx) < EPS && abs(x2[1] - yb) < EPS;
        bool c2 = (cond1 || cond2 || cond3 || cond4 || cond5 ||
                    cond6 || cond7);
        boundaryMask->at((*F)(i, 2)) = c2 ? NodeType::BOUNDARY_FIXED : type;

        cond1 = abs(x3[0] - cx) < EPS && abs(x3[2] - cz) < EPS;
        cond2 = abs(x3[0] - cx) < EPS && abs(x3[2] - zb) < EPS;
        cond3 = abs(x3[0] - xb) < EPS && abs(x3[2] - cz) < EPS;
        cond4 = abs(x3[1] - ya) < EPS && abs(x3[2] - cz) < EPS;
        cond5 = abs(x3[1] - yb) < EPS && abs(x3[2] - cz) < EPS;
        cond6 = abs(x3[0] - cx) < EPS && abs(x3[1] - ya) < EPS;
        cond7 = abs(x3[0] - cx) < EPS && abs(x3[1] - yb) < EPS;
        bool c3 = (cond1 || cond2 || cond3 || cond4 || cond5 ||
                    cond6 || cond7);
        boundaryMask->at((*F)(i, 3)) = c3 ? NodeType::BOUNDARY_FIXED : type;

        assert(c0 + c1 + c2 + c3 <= 1);
      }
    }
  }

  // Sort the removed Id's into reverse order for easy removal
  std::sort(idsToBeRemoved.begin(), idsToBeRemoved.end(), std::greater<int>());

  for (auto id = idsToBeRemoved.begin(); id != idsToBeRemoved.end(); ++id) {
      utils::removeRow(*F, *id);
  }

  Mesh<D> adaptiveMesh(*Vp, *F, *boundaryMask, mon,
      numThreads, rho, tau);

  // Create the solver
  MeshIntegrator<D> solver(dt, adaptiveMesh);

  std::vector<double> Ivals;

  clock_t start = clock();
  double Ih;
  double Ihprev = INFINITY;
  int i;
  for (i = 0; i < nSteps; i++) {
    Ih = solver.step(50, 1e-4);

    Ivals.push_back(Ih);

    double dIdt = (Ih - Ihprev) / dt;
    cout << "d/dt = " << (Ih - Ihprev) / dt << endl;
    cout << "Ih = " << Ih << endl;

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
  string IhDir = outDir + "/Ih.txt";
  adaptiveMesh.outputPoints(pointsOutDir.c_str());
  adaptiveMesh.outputSimplices(triangleOutDir.c_str());
  outputVecToFile(IhDir, Ivals);

  delete Vc;
  delete F;
  delete boundaryMask;
}

template <int D>
void setUpBoxExperiment(string testName, ifstream &inputFile,
    int numThreads, MonitorFunction<D> *mon) {
  int nx, ny, nz, nSteps;
  double dt, tau, rho;
  double xa, xb, ya, yb, za, zb;

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
  nz = 0;
  za = 0;
  zb = 0;
  if (D == 2) {
    inputFile >> nx >> ny;
    nPnts = (nx+1)*(ny+1) + nx*ny;
    inputFile >> xa >> xb >> ya >> yb;

  } else {
    inputFile >> nx >> ny >> nz;
    nPnts = (nx+1)*(ny+1)*(nz+1) + nx*ny*nz;
    inputFile >> xa >> xb >> ya >> yb >> za >> zb;
  }

  // Parameters for the mesh
  std::unordered_map<std::string, double> params;
  params["nx"] = nx;
  params["ny"] = ny;
  params["nz"] = nz;
  params["d"] = D;
  params["rho"] = rho;
  params["tau"] = tau;

  params["xa"] = xa;
  params["xb"] = xb;
  params["ya"] = ya;
  params["yb"] = yb;
  params["za"] = za;
  params["zb"] = zb;

  Eigen::MatrixXd *Vc = nullptr;
  Eigen::MatrixXd *Vp = nullptr;
  Eigen::MatrixXi *F = nullptr;
  vector<NodeType> *boundaryMask = nullptr;

  // Generate the initial mesh
  Vc = new Eigen::MatrixXd(nPnts, D);
  Vp = new Eigen::MatrixXd(nPnts, D);
  if (D == 2) {
    F = new Eigen::MatrixXi(4*nx*ny, D+1);
  } else {
    F = new Eigen::MatrixXi(12*nx*ny*nz, D+1);
  }

  boundaryMask = new vector<NodeType>(Vc->rows());

  utils::generateUniformRectMesh<D>(params, Vp, F, boundaryMask,
      type);

  *Vc = *Vp;
    
  Mesh<D> adaptiveMesh(*Vc, *Vp, *F, *boundaryMask, mon,
      numThreads, rho, tau);

  // Create the solver
  MeshIntegrator<D> solver(dt, adaptiveMesh);

  std::vector<double> Ivals;

  clock_t start = clock();
  double Ih;
  double Ihprev = INFINITY;
  int i;
  for (i = 0; i < nSteps; i++) {
    Ih = solver.step(50, 1e-4);

    Ivals.push_back(Ih);

    double dIdt = (Ih - Ihprev) / dt;
    cout << "d/dt = " << (Ih - Ihprev) / dt << endl;
    cout << "Ih = " << Ih << endl;

    if (i != 0 && (abs(dIdt) < 1e-4 || dIdt > 0)) {
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
  string Ihdir = outDir + "/Ih.txt";
  outputVecToFile(Ihdir, Ivals);

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

  int D;
  inputFile >> D;

  // Specify the monitor function
  auto *M0 = new MEx0<2>();
  auto *M1 = new MEx1<2>();
  auto *M2 = new MEx2<2>();
  auto *M3 = new MEx3<2>();
  auto *M03D = new MEx0<3>();
  auto *M13D = new MEx13D<3>();
  auto *M23D = new MEx23D<3>();
  auto *M33D = new MEx33D<3>();

  vector<MonitorFunction<2>*> Mvals;
  vector<MonitorFunction<3>*> Mvals3D;
  Mvals.push_back(M0);
  Mvals.push_back(M1);
  Mvals.push_back(M2);
  Mvals.push_back(M3);

  Mvals3D.push_back(M03D);
  Mvals3D.push_back(M13D);
  Mvals3D.push_back(M23D);
  Mvals3D.push_back(M33D);

  // Get monitor id
  int monType;
  inputFile >> monType;

  if (testType.compare("SquareGrid") == 0) {
    if (D == 2) {
      setUpBoxExperiment<2>(inFileName, inputFile, numThreads, Mvals.at(monType));
    } else {
      setUpBoxExperiment<3>(inFileName, inputFile, numThreads, Mvals3D.at(monType));
    }
  } else if (testType.compare("LevelSet") == 0) {
    if (D == 2) {
      setUpLevelSetExperiment<2>(inFileName, inputFile, numThreads, Mvals.at(monType));
    } else {
      setUpLevelSetExperiment<3>(inFileName, inputFile, numThreads, Mvals3D.at(monType));
    }
  } else if (testType.compare("Shoulder") == 0) {
    if (D == 2) {
      setUpShoulderExperiment<2>(inFileName, inputFile, numThreads, Mvals.at(monType));
    } else {
      setUpShoulderExperiment<3>(inFileName, inputFile, numThreads, Mvals3D.at(monType));
    }
  } else {
    assert(false);
  }

  inputFile.close();
}
