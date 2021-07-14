#include "AdaptationFunctional.h"
#include <Eigen/Dense>
#include "MonitorFunction.h"
#include <vector>
#include <iostream>
#include <Eigen/StdVector>

using namespace std;

template <int D>
AdaptationFunctional<D>::AdaptationFunctional(const AdaptationFunctional &obj) {
    Vc = obj.Vc;
    Vp = obj.Vp;
    DXpU = obj.DXpU;

    this->M = obj.M;
    this->boundaryNode = obj.boundaryNode;
    cout << "In CTOR" << endl;
    assert(false);

    // assert(false);
    mPre = new vector<Eigen::Matrix<double, D, D>>(D+1);
    for (int i = 0; i < D+1; i++) {
        mPre->at(i) = obj.mPre->at(i);
    }
    // for (int i = 0; i < D+1; i++) {
    //     // Eigen::Matrix<double, D, D> temp = Eigen::Matrix<double, D, D>::Zero();
    //     // mPre->push_back(temp);
    // }
    // assert(false);
}

/**
 * Bring the specified id to the front of the vector via a swap
*/
template <int D>
void swap(int id, Eigen::Vector<int, D+1> &ids) {
    for (int i = 0; i < ids.size(); i++) {
        if (ids[i] == id) {
            ids[i] = ids[0];
            ids[0] = id;
            break;
        }
    }
}

template <int D>
AdaptationFunctional<D>::AdaptationFunctional( Eigen::MatrixXd &Vc,
        Eigen::MatrixXd &Vp, Eigen::MatrixXi &F,
        Eigen::VectorXd &DXpU, MonitorFunction<D> *m, double w) {
    
    if (D == -1) {
        cout << "Must specify DIM >= 1!" << endl;
        assert(D == -1);
    }
    
    this->w = w;
    this->boundaryNode = boundaryNode;

    this->DXpU = &DXpU;

    this->Vc = &Vc;
    this->Vp = &Vp;
    this->DXpU = &DXpU;

    this->M = m;

    this->mPre = new vector<Eigen::Matrix<double, D, D>>(D+1);
}


/**
 * Method for computing the block gradient of a simplex for the objective function
*/
template <int D>
double AdaptationFunctional<D>::blockGrad(int zId, Eigen::Vector<double, D*(D+1)> &z,
            Eigen::Vector<double, D*(D+1)> &xi,
            Eigen::Vector<double, D*(D+1)> &grad,
            MeshInterpolator<D> &interp) {
    double Ih = 0.0;;
    double detFJ;
    Eigen::Vector<double,D> gradSimplex;

    Eigen::Matrix<double,D,D> E;
    Eigen::Matrix<double,D,D> FJ;
    Eigen::Matrix<double,D,D> Ehat;
    Eigen::Matrix<double,D,D> M;
    Eigen::Vector<double,D> xK;
    Eigen::Matrix<double,D,D> vLoc;
    Eigen::Vector<double,D> xTemp;
    Eigen::Matrix<double,D,D> Einv;
    Eigen::Matrix<double,D,D> dGdJ;
    Eigen::Matrix<double,D,D> Mt0;
    Eigen::Matrix<double,D,D> Mtn;
    Eigen::Vector<double,D> dGdX;
    Eigen::Matrix<double,D,D> dGdM;
    Eigen::Vector<double,D> basisComb;
    Eigen::Vector<int, D+1> ids;

    double dFact;
    if (D == 2) {
        dFact = 2.0;
    } else if (D == 3) {
        dFact = 6.0;
    }
    double absK;

    // Build the centroid
    xK = z.segment(0, D);
    for (int n = 1; n <= D; n++) {
        xK += z.segment(D*n, D);
    }
    xK /= ((double) D + 1.0);

    // Interpolate the monitor function
    // cout << "zId = " << zId << endl;
    Eigen::Vector<double, D+1> bCoords;
    // int sId = interp.evalWithKnn(xK, bCoords);
    int sId = zId;
    interp.evalMonitorOnSimplex(sId, xK, bCoords, M);
    // interp.evalMonitorOnSimplex(zId, xK, M);
    Eigen::Matrix<double, D, D> Minv(M.inverse());
    // Eigen::Matrix<double, D, D> Minv;

    for (int i = 0; i < D+1; i++) {
        Eigen::Matrix<double, D, D> mTemp;
        xTemp = z.segment(i*D, D);
        // sId = interp.evalWithKnn(xTemp, bCoords);
        interp.evalMonitorOnSimplex(sId, xTemp, bCoords, mTemp);
        mPre->at(i) = mTemp;
    }

    double G, dGddet;

    for (int i = 0; i < D+1; i++) {
        // Compute the edge matrix
        int j = 0;
        for (int n = (i+1)%(D+1); n != i; n = (n+1)%(D+1)) {
            E.col(j)    = z.segment(D*n, D)  - z.segment(D*i, D);
            Ehat.col(j) = xi.segment(D*n, D) - xi.segment(D*i, D);
            j++;
        }

        // The element volume
        double Edet = E.determinant();

        Einv = E.inverse();

        if (i == 0) {
            // Approximate Jacobian
            FJ = Ehat * Einv;
            detFJ = FJ.determinant();
        }

        absK = abs(Edet/dFact);

        if (i == 0) {
            G = this->G(FJ, detFJ, Minv, xK);
            this->dGdJ(FJ, detFJ, Minv, xK, dGdJ);
            dGddet = this->dGddet(FJ, detFJ, Minv, xK);
            this->dGdX(FJ, detFJ, Minv, xK, dGdX);
            this->dGdM(FJ, detFJ, Minv, xK, dGdM);
        }
        
	
        basisComb.setZero();
        j = 0;
        for (int n = (i+1)%(D+1); n != i; n = (n+1)%(D+1)) {
            basisComb += Einv.row(j) * ((dGdM * (mPre->at(n) - mPre->at(i))).trace());
            j++;
        }

        vLoc = -G*Einv + Einv*dGdJ*Ehat*Einv + dGddet*detFJ*Einv;
        for (int n = 0; n < D; n++) {
            vLoc.row(n) -= (basisComb + dGdX)/((double) D + 1.0);
        }

        // Compute the gradient for current simplex
        gradSimplex.setZero();

        for (int n = 0; n < D; n++) {
            gradSimplex += vLoc.row(n);
        }
        gradSimplex += basisComb + dGdX;

        gradSimplex *= absK;

        // Compute the gradient
        for (int l = 0; l < D; l++) {
            grad(D*i + l) = gradSimplex(l);
        }

        // Update the energy
        if (i == 0) {
            Ih += absK * G;
        }
    }

    // Now add the constraint regularization
    // Ih += 0.5*w*w*( (*DXpU).segment(D*(D+1)*zId, D*(D+1)) - z ).squaredNorm();

    grad += -w*w*(*DXpU).segment(D*(D+1)*zId, D*(D+1)) + w*w*z;

    return Ih;
}



/**
 * Method for computing the block gradient of a simplex for the objective function
*/
template <int D>
double AdaptationFunctional<D>::blockGradC(int zId, Eigen::Vector<double, D*(D+1)> &z,
            Eigen::Vector<double, D*(D+1)> &xi,
            Eigen::Vector<double, D*(D+1)> &grad,
            MeshInterpolator<D> &interp) {
    double Ih = 0.0;;
    double detFJ;
    double gradSimplex[D] = {0};  
    Eigen::Matrix<double,D,D> E;
    // Eigen::Matrix<double,D,D> FJ;
    Eigen::Matrix<double,D,D> Ehat;
    Eigen::Matrix<double,D,D> M;
    Eigen::Vector<double,D> xK;
    Eigen::Matrix<double,D,D> vLoc;
    Eigen::Vector<double,D> xTemp;
    // Eigen::Matrix<double,D,D> Einv;
    Eigen::Matrix<double,D,D> dGdJ;
    double dGdJ1[D][D] = {0};
    double dGdX1[D] = {0};
    double dGdM1[D][D] = {0};
    double basisComb[D] = {0};


    double dFact;
    if (D == 2) {
        dFact = 2.0;
    } else if (D == 3) {
        dFact = 6.0;
    }
    double absK;

    // Build the centroid
    xK = z.segment(0, D);
    for (int n = 1; n <= D; n++) {
        xK += z.segment(D*n, D);
    }
    xK /= ((double) D + 1.0);

    // Interpolate the monitor function
    Eigen::Vector<double, D+1> bCoords;
    int sId = zId;
    interp.evalMonitorOnSimplex(sId, xK, bCoords, M);
    double Minv[D][D] = {0};
    Boost_Inverse(M, Minv);

    for (int i = 0; i < D+1; i++) {
        Eigen::Matrix<double, D, D> mTemp;
        xTemp = z.segment(i*D, D);
        // sId = interp.evalWithKnn(xTemp, bCoords);
        interp.evalMonitorOnSimplex(sId, xTemp, bCoords, mTemp);
        mPre->at(i) = mTemp;
    }

    double G, dGddet;

    double EinvC[D][D];
    // double EhatC[D][D];
    double FJC[D][D];
    for (int i = 0; i < D+1; i++) {
        // Compute the edge matrix
        int j = 0;
        for (int n = (i+1)%(D+1); n != i; n = (n+1)%(D+1)) {
            E.col(j)    = z.segment(D*n, D)  - z.segment(D*i, D);
            Ehat.col(j) = xi.segment(D*n, D) - xi.segment(D*i, D);
            j++;
        }

        // The element volume
        double Edet = E.determinant();

        Boost_Inverse(E, EinvC);
        // Boost_cast_EtoS(Ehat, EhatC);

        // Einv = E.inverse();

        if (i == 0) {
            // Approximate Jacobian
            // FJ = Ehat * Einv;
            Boost_Multi(Ehat, EinvC, FJC);
            // detFJ = FJ.determinant();
            detFJ = Boost_determinant(FJC);
        }

        absK = abs(Edet/dFact);

        if (i == 0) {
            // Evaluating the mesh functional
            double FJ1[D][D] = {0};
            double xK1[D] = {0};
            for(int n = 0; n < D; n++){
                 xK1[n] = xK(n);
                 for(int l = 0; l < D; l++){
                     FJ1[n][l] = FJC[l][n];
                 }
            }
            G = this->GC(FJ1,detFJ, Minv, xK1);
            this->dGdJC(FJ1,detFJ,Minv,xK1,dGdJ1);
            dGddet = this->dGddetC(FJ1, detFJ, Minv, xK1);
            this->dGdXC(FJ1, detFJ, Minv, xK1, dGdX1);
            this->dGdMC(FJ1, detFJ, Minv, xK1, dGdM1);
        }
        
        // basisComb.setZero(); // Alrady set to zero
        j = 0;
        // Eigen::Vector<double, D> bTemp;
        double basisCombC[D] = {0};
        Eigen::Matrix<double, D,D> mpretemp;
        // bTemp.setZero();
        for (int n = (i+1)%(D+1); n != i; n = (n+1)%(D+1)) {
            // Computing bTemp += Einv.row(j) * ((dGdM * (mPre->at(n) - mPre->at(i))).trace());
            mpretemp = (mPre->at(n) - mPre->at(i));
            double dgdm_mpre[D][D] = {0};
            Boost_Multi(dGdM1, mpretemp, dgdm_mpre);
            double tr = Boost_Trace(dgdm_mpre);

            for (int k = 0; k < D; k++) {
                basisCombC[k] += EinvC[j][k]*tr;
            }
            j++;
        }

        // Setting vLoc = -G*Einv + Einv*dGdJ*Ehat*Einv + dGddet*detFJ*Einv;
        // Computing first part: vLoc = -G*Einv;
        for (int n = 0; n < D; n++) {
            for (int m = 0; m < D; m++) {
                vLoc(n, m) = -G*EinvC[n][m];
            }
        }

        // Computing: Einv*dGdJ
        double tmp[D][D] = {0};
        Boost_Multi(EinvC, dGdJ1, tmp);

        // Computing: Einv*dGdJ*Ehat*Einv + dGddet*detFJ*Einv;
        double t1[D][D] = {0};
        double t2[D][D] = {0};
        Boost_Multi(tmp, Ehat, t1);
        Boost_Multi(t1, EinvC, t2);
        for (int n = 0; n < D; n++) {
            for (int m = 0; m < D; m++) {
                vLoc(n, m) += t2[n][m] + dGddet * detFJ * EinvC[n][m];
            }
        }

        for (int n = 0; n < D; n++) {
            for (int l = 0; l < D; l++){
                vLoc(n,l) -= (basisCombC[l] + dGdX1[l])/((double) D + 1.0);//0.748216
            }
        }

        // Compute the gradient for current simplex
        Boost_row(vLoc, gradSimplex);

        for (int n = 0; n < D; n++){
            gradSimplex[n] += basisCombC[n] + dGdX1[n];
        }

        for (int n = 0; n < D; n++){
            gradSimplex[n] *= absK;
        }
        // Compute the gradient
        for (int l = 0; l < D; l++) {
            grad(D*i + l) = gradSimplex[l];
        }

        // Update the energy
        if (i == 0) {
            Ih += absK * G;
        }
    }

    // Now add the constraint regularization
    // Ih += 0.5*w*w*( (*DXpU).segment(D*(D+1)*zId, D*(D+1)) - z ).squaredNorm();

    grad += -w*w*(*DXpU).segment(D*(D+1)*zId, D*(D+1)) + w*w*z;

    return Ih;
}
template <int D>
AdaptationFunctional<D>::~AdaptationFunctional() {
    delete mPre;
}

template class AdaptationFunctional<2>;
template class AdaptationFunctional<3>;
