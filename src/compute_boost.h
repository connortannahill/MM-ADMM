#ifndef COMPUTE_BOOST_H
#define COMPUTE_BOOST_H

#include <iostream>
#include <Eigen/Dense>

// Compute Determinant of Eigen matrix. 
// Not efficient.
template<int D=-1>
inline double Boost_determinant(Eigen::Matrix<double,D,D> &M){
    return M.determinant();
}

// Compute Determinant of static matrix.
template<int D=-1>
inline double Boost_determinant(double (&M)[D][D]){
    if(D == 2){
        return M[0][0]*M[1][1] - M[1][0]*M[0][1];
    } else {
        double tmp1, tmp2, tmp3;
        tmp1 = M[1][1]*M[2][2] - M[1][2]*M[2][1];
        tmp2 = M[1][0]*M[2][2] - M[1][2]*M[2][0];
        tmp3 = M[1][0]*M[2][1] - M[1][1]*M[2][0];
        return M[0][0]*tmp1 - M[0][1]*tmp2 + M[0][2]*tmp3;
    }
}

// Compute Inverse of an Eigen matrix and store to static array.
template<int D=-1>
inline void Boost_Inverse(Eigen::Matrix<double,D,D> &M, double (&Minv)[D][D]){
    double detM = 0;
    detM = M.determinant();
    double tmp[D][D] = {0};
    Eigen::Map<Eigen::Matrix<double,D,D>>(&tmp[0][0], D, D) = M;
    if(D == 2){
        Minv[0][0] = tmp[1][1] / detM;
        Minv[0][1] = -1*tmp[1][0] / detM;
        Minv[1][0] = -1*tmp[0][1] / detM;
        Minv[1][1] = tmp[0][0] / detM;
    } else {
        for(int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
              //Minv[i][j] = ((tmp[(j+1)%3][(i+1)%3] * tmp[(j+2)%3][(i+2)%3]) - (tmp[(j+1)%3][(i+2)%3] * tmp[(j+2)%3][(i+1)%3])) / detM;
              Minv[i][j] = ((tmp[(i+1)%3][(j+1)%3] * tmp[(i+2)%3][(j+2)%3]) - (tmp[(i+1)%3][(j+2)%3] * tmp[(i+2)%3][(j+1)%3])) / detM;
            }
        }
    }
    return;
}

// Multiplication function. result = constant * (A x B). B_transpose determine if A x B^T.  
template<int D = -1>
inline void Boost_Multi(double (&A)[D][D], double (&B)[D][D], double (&result)[D][D],bool B_transpose=false, double constant=1.0){
    if(!B_transpose){
        for(int i = 0; i < D; i++){
            for (int j = 0; j < D; j++){
                double xg = 0;
                for(int o = 0; o < D; o++){
                    xg += A[i][o] * B[o][j];
                }
                result[i][j] = constant*xg;
            }
        }
    } else {
        for(int i = 0; i < D; i++){
            for (int j = 0; j < D; j++){
                double xg = 0;
                for(int o = 0; o < D; o++){
                    xg += A[i][o] * B[j][o];
                }
                result[i][j] = constant*xg;
            }
        }
    }
    return;
}

template<int D = -1>
inline void Boost_cast_EtoS(Eigen::Matrix<double, D, D> &M, double (&Mcast)[D][D]){
    Eigen::Map<Eigen::Matrix<double, D, D>>(&Mcast[0][0], D, D) = M;
    for(int i = 0; i < D; ++i){
        for(int j = i; j < D; ++j){
            double tmp = Mcast[i][j];
            Mcast[i][j] = Mcast[j][i];
            Mcast[j][i] = tmp;
        }
    }
    return;
}

template<int D = -1>
inline void Boost_Multi(Eigen::Matrix<double, D,D> &M, double (&B)[D][D], double (&result)[D][D],bool B_transpose=false, double constant=1.0){
    double A[D][D] = {0};
    Eigen::Map<Eigen::Matrix<double, D,D>>(&A[0][0],D,D) = M;
    // Boost_cast_EtoS(M, A);
    if(!B_transpose){
        for(int i = 0; i < D; i++){
            for (int j = 0; j < D; j++){
                double xg = 0;
                for(int o = 0; o < D; o++){
                    xg += A[o][i] * B[o][j];
                }
                result[i][j] = constant*xg;
            }
        }
    } else {
        for(int i = 0; i < D; i++){
            for (int j = 0; j < D; j++){
                double xg = 0;
                for(int o = 0; o < D; o++){
                    xg += A[o][i] * B[j][o];
                }
                result[i][j] = constant*xg;
            }
        }
    }
    return;
}

template<int D = -1>
inline void Boost_Multi(double (&A)[D][D], Eigen::Matrix<double, D, D> &M, double (&result)[D][D],bool B_transpose=false, double constant=1.0){
    double B[D][D] = {0};
    Eigen::Map<Eigen::Matrix<double, D,D>>(&B[0][0],D,D) = M;
    if(!B_transpose){
        for(int i = 0; i < D; i++){
            for (int j = 0; j < D; j++){
                double xg = 0;
                for(int o = 0; o < D; o++){
                    xg += A[i][o] * B[j][o];
                }
                result[i][j] = constant*xg;
            }
        }
    } else {
        for(int i = 0; i < D; i++){
            for (int j = 0; j < D; j++){
                double xg = 0;
                for(int o = 0; o < D; o++){
                    xg += A[i][o] * B[o][j];
                }
                result[i][j] = constant*xg;
            }
        }
    }
    return;
}

template<int D = -1>
inline double Boost_Trace(double (&M)[D][D]){
    double ans = 0;
    for(int i = 0; i < D; ++i){
        ans += M[i][i];
    }
    return ans;
}

template<int D = -1>
inline void Boost_row(Eigen::Matrix<double, D, D> &M, double (&result)[D]){
    double tmp[D][D] = {0};
    Eigen::Map<Eigen::Matrix<double, D, D>>(&tmp[0][0],D,D) = M;
    for (int n = 0; n < D; n++) {
        for (int l = 0; l < D; l++){
            result[l] += tmp[l][n];
        }
    }
}

template<int D = -1>
inline void Boost_row(double (&A)[D][D], double (&result)[D]){
    //double tmp[D][D] = {0};
    //Eigen::Map<Eigen::Matrix<double, D, D>>(&tmp[0][0],D,D) = M;
    for (int n = 0; n < D; n++) {
        for (int l = 0; l < D; l++){
            result[l] += A[n][l];
        }
    }
}


template<int D = -1>
inline void Boost_col(Eigen::Matrix<double, D, D> &M, double (&result)[D]){
    double tmp[D][D] = {0};
    Eigen::Map<Eigen::Matrix<double, D, D>>(&tmp[0][0],D,D) = M;
    for (int n = 0; n < D; n++) {
        for (int l = 0; l < D; l++){
            result[l] += tmp[n][l];
        }
    }
}

template<int D = -1>
inline void Boost_Add(double (&A)[D][D], double (&B)[D][D], double (&result)[D][D],double scalar=1.0){
    for(int i = 0; i < D; ++i){
        for(int j = 0; j < D; ++j){
            result[i][j] = scalar*(A[i][j] + B[i][j]);
        }
    }
}


template<int D = -1>
inline void Boost_Minus(double (&A)[D][D], double (&B)[D][D], double (&result)[D][D],double scalar=1.0){
    for(int i = 0; i < D; ++i){
        for(int j = 0; j < D; ++j){
            result[i][j] = scalar*(A[i][j] - B[i][j]);
        }
    }
}

#endif