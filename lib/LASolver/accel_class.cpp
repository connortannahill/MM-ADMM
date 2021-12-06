#include "def_compiler.h"
#include "Standard.h"
#include "SparseItUtil.h"
#include "ILU_class.h"
#include "accel_class.h"

using namespace SparseItObj;
using namespace std;

namespace SparseItObj { // start namespace

scaler_orthomin::scaler_orthomin(int n_in, // number of unknowns
                                 double* res_in, // initial residual
                                 int north_in // number of orthog vectors
                                )
{
    n = n_in;
    north = north_in;

    int i;
    rmsi = 0.0;

    for (i=0; i<n; i++) {
      rmsi += res_in[i]*res_in[i];
    }// end initialzation
    rmsi = sqrt( rmsi );

    vk = NULL;
    q = NULL;
    Aq = NULL;
    Avk = NULL;
    AqAq = NULL;
    a_mult = NULL;

    vk = new double [n];
    assert (vk != NULL);

    Avk = new double [n];
    assert( Avk != NULL);

    double* dptr1 = NULL;
    double* dptr2 = NULL;

    try{// start try block

      dptr1 = new double [n*north];
      assert( dptr1 != NULL);
      dptr2 = new double [n*north];
      assert( dptr2 != NULL);

      Aq = new double* [north];
      assert( Aq != NULL);
      q = new double* [north];
      assert( q != NULL);

      AqAq = new double [north];
      assert( AqAq != NULL);
      a_mult = new double[north];
      assert( a_mult!= NULL);

    }// end try block
    catch(std::bad_alloc) {
      delete [] dptr1; dptr1 = NULL;
      delete [] dptr2; dptr2 = NULL;
      throw;
    }
    

    {// start block

      int count = 0;

      for(i=0; i< north; i++){ // start loop
        Aq[i]  = &dptr1[count];
        q[i] =   &dptr2[ count];
        count += n;
      }// end loop

   } // end block

   iter_count = 0; // no previous orthog vectors

}// end scaler_orthomin constructor

scaler_orthomin::~scaler_orthomin(void)
{
   delete [] vk; vk = NULL;
   delete [] Avk; Avk = NULL;
     
   if (Aq != NULL) {
     delete [] Aq[0];
   }

   delete [] Aq; Aq = NULL;

   if (q != NULL) {
     delete [] q[0];
   }
   delete [] q; q = NULL;
   delete [] AqAq; AqAq = NULL;
   delete [] a_mult; a_mult = NULL;

}// end scaler_orthomin destructor

void scaler_orthomin::acc_scaler(int ia[], int ja[], double A[], 
                  int iter,  // iteration count
                  double sol[],
                          // on entry, x^0 s.t. res^0 =
                          //              b - A x^0
                          // on exit, sol = A^{-1} res^0 + x^0
                  double res[],
                  double ctol, // residual reduction
                  int *iconv,  //  = 0 converged
                               // !=0 not converged
                  const double toler[], // change convergence criteria
                  double *rms,  // output value of rms error
                  scaler_ILU &ilu)
{
  const double tiny = 1.e-300;
     
  if (ilu.check_factor() == 1) {
    throw General_Exception("Error numeric ILU not done orthomin\n");
  }

  ilu.solve(vk, res); // LU vk = res

  matmult( vk, Avk, n, A, ia, ja); 

  X1EqualsX2( q[iter_count], vk, n);
  X1EqualsX2( Aq[iter_count], Avk, n);

  if (iter_count > 0) { // previous vectors

    int i;
    for (i=0; i<= iter_count-1; i++) {// loop over previous
      a_mult[i] = - dot(Avk, Aq[i],n )/( AqAq[i] + tiny);
    }// end loop over previous

    for (i=0; i<= iter_count-1; i++) {// loop over previous

      X1EqualsX1PlusAlphaTimesX2(q[iter_count], a_mult[i], q[i], n);
      X1EqualsX1PlusAlphaTimesX2(Aq[iter_count], a_mult[i], Aq[i], n);

    }// end loop over previous
         
  }// end previous vectors

  // minimization
  AqAq[iter_count] = dot( Aq[iter_count], Aq[iter_count], n);

  double ohm = dot( Aq[iter_count], res, n)/( AqAq[iter_count] + tiny);

  {// start update block
    int i;
    *rms = 0.0;
    *iconv = 0;

    double step;
    double temp;

    for (i=0; i<n; i++) { // loop over unknowns
      step = ohm*q[iter_count][i];
      sol[i] += step;

      res[i] -=  ohm*Aq[iter_count][i];
  
      // convergence testing
      temp = max( fabs(step), fabs(vk[i]) );
      temp = max( temp   , fabs( q[iter_count][i] ) );
      if (temp > fabs(toler[i])) {
        (*iconv)++;
      }
      *rms = *rms + res[i]*res[i];
    }// end loop over unknowns

    *rms = sqrt(*rms);
    if ((*rms/rmsi) <  ctol) {
      *iconv = 0;
    }

  }// end update block
     
  iter_count++;

  if (iter_count > north-1) { // restart
    iter_count = 0;
  }// end restart

  return;

}// end scaler_orthomin::acc_scaler

scaler_conj::scaler_conj(int nin, double *res_in) : base_accel()
{
  n = nin;
  int i;

  avk = NULL;
  aq = NULL;
  q = NULL;
  v = NULL;

  avk = new double [n];
  assert( avk != NULL);
  aq = new double [n];
  assert( aq != NULL);
  q = new double [n];
  assert( q != NULL);
  v = new double [n];
  assert( v != NULL);

  rmsi = 0.0;
  for (i = 0; i <= n - 1; i++) {
    rmsi = rmsi + res_in[i]*res_in[i];
    q[i] = 0.0;
  }
  rmsi = sqrt(rmsi);

} // end scaler_conj( int nin, double *res_in)

scaler_conj::~scaler_conj(void)
{
  delete []  aq ;
  delete [] avk ;
  delete [] q ;
  delete [] v ;

} // end ~scaler_conj(void)

void scaler_cgstab::init( double *res_in)
{
  int i;
  alpha = 1.0;
  rholst = 1.0;
  omega = 1.0;

  res0 = NULL; pvec = NULL; vbar = NULL;
  avbar = NULL; svec = NULL; zvec = NULL;
  tvec = NULL; workx = NULL;

  res0  = new double [n];
  assert( res0 != NULL);
  pvec  = new double [n];
  assert( pvec != NULL);
  vbar  =  new double [n];
  assert( vbar != NULL);
  avbar  = new double [n];
  assert( avbar != NULL);
  svec  =  new double [n];
  assert( svec != NULL);
  zvec  = new double [n];
  assert( zvec != NULL);
  tvec  = new double [n];
  assert( tvec != NULL);
  workx  = new double [n];
  assert( workx != NULL);

  rmsi = 0.0;
  for (i=0;i <= n-1; i++) {
    pvec[i] = 0.0;
    avbar[i] = 0.0;
    res0[i] = res_in[i];
    rmsi += res_in[i] * res_in[i];
  }
  rmsi = sqrt(rmsi);
}

scaler_cgstab::~scaler_cgstab(void)
{
  delete [] res0;
  delete [] pvec;
  delete [] vbar;
  delete [] avbar;
  delete [] svec;
  delete [] zvec;
  delete [] tvec ;
  delete [] workx ;
}

void scaler_cgstab::acc_scaler(int ia[], int ja[], double a[], int iter, 
                               double sol[], 
                                       // on entry, x^0 s.t. res^0 = 
                                       //              b - A x^0
                                       // on exit, sol = A^{-1} res^0 + x^0
                               double res[], double ctol, int *iconv, 
                               const double toler[], double *rms, 
                               scaler_ILU &ilu)
/*

  solve a * sol = res
    use cgstab acceleration

    input:  
       ia, ja
       a[]
       iter    iteration number
       res[]      init right hand side
       ctol   resid reduction tol
       itmax      max number of iterations
       ilu    ref to ILU class
       n      class variable

       *res0, *pvec, *vbar, *avbar, *svec, *zvec, *tvec, *workx, rmsi;
       omega, rholst, alpha;

       class wk variables, initialized in 
       constructor, see "accel_class.h"

   output:
    *iconv   =0  converged
            !=0  not converged

    sol[]     solution
    *rms      rms resid

*/
{
  /*   
    cgstab acceleration
  */
  double beta, rho, step, step2, step1,  tiny;
  int i;

  if (ilu.check_factor() == 1) {
    throw General_Exception("Error numeric ILU not done in cgstab\n");
  }

  tiny = 1.e-300;

/*
 rho = (res0, res)
 beta = (rho/rholst)*(alpha/omega)
 rholst = rho
 pvec = res + beta*(pvec - omega*Avbar)
*/

  rho = dot(res0 , res , n);
  beta = rho / (rholst + tiny);
  beta *= (alpha/(omega + tiny));
  rholst = rho;
  for (i=0; i<= n-1; i++) {
    pvec[i] = res[i] + beta * (pvec[i]-omega*avbar[i]);
  }

/*
 solve (LU) vbar = pvec
 Avbar = A*vbar
 alpha = rholst/(res0,Avbar)
*/

  ilu.solve(vbar, pvec);
  matmult(vbar, avbar, n, a, ia, ja);
  alpha = rho / dot(res0 , avbar, n);

/*
 svec = res - alpha*Avbar
 solve  (LU) zvec = svec
 Azvec = A*zvec = tvec
 omega = (A zvec,svec)/(A zvec,A zvec) in CG-STAB, but
 omega = (tvec,svec)/(tvec,tvec) in CG-STAB-P
*/

  for (i=0; i<=n-1; i++) {
    svec[i] = res[i] - alpha * avbar[i];
  }
  ilu.solve(zvec, svec);
  matmult(zvec, tvec, n, a, ia, ja);
  omega = dot(tvec, svec, n) / (dot(tvec, tvec, n) + tiny);

/*
 sol = sol + alpha*vbar + omega*zvec
 res = svec - omega*tvec
*/

  *rms = 0.0;
  *iconv = 0;
  for (i=0; i<= n-1; i++) {
    step1 = alpha*vbar[i];
    step2 = omega*zvec[i];
    step = step1 + step2;
    if (fabs(step) > fabs(toler[i])) {
      (*iconv)++;
    }
      sol[i] = sol[i] + step;
      res[i] = svec[i] - omega*tvec[i];
      *rms = *rms + res[i]*res[i];
  }

  *rms = sqrt( *rms );
  if ((*rms/rmsi) <  ctol) {
    *iconv = 0;
  }

  return;

/*
  end of cgstab
*/
}


void scaler_conj::acc_scaler(int ia[], int ja[], double a[], int iter, 
double x[], double res[], double ctol, int *iconv, const double toler[],
double *rms, scaler_ILU &ilu)
 
/*     conj gradient acceleration     */

/*   input:
       ctol      conver tol (rms/rms_0 < ctol)
       a[], ia, ja   orig matrix
       res[]     init residual
       iter      iteration number
       x[]       current solutin guess
       ilu ref to ILU class
       n     class variable
      toler   chages conv tol vector

    *q, *aq, *avk, *v;
     rvkm, rmsi;
           class variables, initialized in constructor
           see "accel_class.h"

    output:
      res[]   updated resid
      *iconv    =0 converged
                !=0 not converged
      x[]     updated soln
      *rms    current rms
*/
{
  int i;
  double rvk, aconj, ohm;

  if (ilu.check_factor() == 1) {
    throw General_Exception("Error numeric ILU not done in conj\n");
  }

  for (i = 0; i <= n-1; i++) {
    avk[i] = 0.0;
  }

  /* solver LU v = res */
  ilu.solve(v, res);

  /* old order */
  /*
    avk = A * v
  */
  matmult(v , avk, n, a, ia, ja);

  /*  
    rvk = (res, v)
  */
  rvk = dot(res, v, n);
 
  if (iter == 1) {
    X1EqualsX2(q, v, n);
    X1EqualsX2(aq, avk, n);
  }
  else {
    aconj = rvk / rvkm;
    for (i = 0; i<= n-1; i++) {
      q[i] = v[i] + aconj * q[i];
      aq[i] = avk[i] + aconj * aq[i];
    }
  }
 
  ohm = rvk/ dot( q, aq, n);
  rvkm = rvk;
 
  *rms = 0.0;
  *iconv = 0;
  for (i = 0; i<= n -1; i++) {
    x[i] =  x[i] + ohm * q[i];
    res[i] = res[i] - ohm * aq[i];
    *rms = *rms + res[i]*res[i];
    if (fabs(ohm*q[i]) > fabs(toler[i]) ||
        fabs(q[i]) > fabs(toler[i]) ||
        fabs(v[i]) > fabs(toler[i]) )
    {
      (*iconv)++;
    }
  }

  *rms = sqrt( *rms );
  if ((*rms/rmsi) <  ctol) {
    *iconv = 0;
  }

  return;
}


//   non member functions
 
double dot(double x1[], double x2[], int n)
/*
     dot product of two vectors
 
input:
    x1[]
    x2[]
    n    length of x1, xe
 
output:
 
    ( x1, x2 )
*/
{
  double temp;
  int i ;
  temp = 0.0;
  for (i=0; i < n ; i++) {
    temp +=  x1[i] * x2[i];
  }

  return temp;

}

void matmult(double xin[], double xout[], int n, double a[], int ia[], 
             int ja[])
/*
  matrix vec multiply
 
input:
 
    a[], ia, ja    matrix in ysmp format
 
    xin[]          input vector
 
output:
 
    xout = A * xin
*/
{
  int i,ii;
  int upper;
  int lower;

  for (i=0; i<n; i++) {
    xout[i] = 0.0;
    lower = ia[i];
    upper = ia[i+1];
    for (ii =lower ; ii <upper; ii++) {
      xout[i] += a[ii] * xin[ ja[ii] ];
    }
  }
}

void X1EqualsX2(double* x1, double* x2, int n)
{
  int i;

  for (i=0; i<n; i++) {
    x1[i] = x2[i];
  }

  return;

}//end X1EqualsX2


void X1EqualsX1PlusAlphaTimesX2( double* x1, double alpha, double* x2, int n)
{
  int i;

  for (i=0; i<n; i++) {
    x1[i] += alpha * x2[i];
  }

}// end X1PlusAlphaTimesX2


} // end namespace SparseItObj
