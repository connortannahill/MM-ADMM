#ifndef ACCEL_CLASS_INC
#define ACCEL_CLASS_INC


namespace SparseItObj { // start namespace SparseItObj

class base_accel
{
  public:

    // default constructor
    base_accel(void){ };
 
    // virtual destructor
    virtual ~base_accel(void){}
 
    virtual void acc_scaler(int ia[], int ja[], double a[], int iter, 
                            double sol[], double res[], double ctol, 
                            int *iconv, const double toler[], double *rms, 
                            scaler_ILU &ilu) = 0; 
}; // end base_accel


class scaler_orthomin : public base_accel
{
  private:

    int n; // number of unknowns
    int north; // max number of orthognel vectors

    double *vk, **q, **Aq, *Avk;
    double *AqAq, *a_mult;
    double rmsi;

    int iter_count; // count of itns since last restart

    // constructors

  public:

    scaler_orthomin(void) {
      cout << "default constructor for orthomin:error \n";
      exit(1);
    }

    scaler_orthomin(int n_in, // number of unknowns
                    double* res_in, // initial residual
                    int north_in // number of orthog vectors
                   );

    // destructor
    virtual ~scaler_orthomin(void);

    virtual void acc_scaler(int ia[], int ja[], double a[], int iter, 
                            double sol[], double res[], double ctol, 
                            int *iconv, const double toler[], double *rms, 
                            scaler_ILU &ilu);
                            // on entry, x^0 s.t. res^0 = b - A x^0
                            // on exit, sol = A^{-1} res^0 + x^0
};// end scaler_orthomin


class scaler_cgstab : public base_accel
{
  private:

    int n;
    double *res0, *pvec, *vbar, *avbar, *svec, *zvec, *tvec, *workx, rmsi;
    double omega, rholst, alpha;

    // constructors

  public:

    scaler_cgstab(void) {
      cout<<"default constructor for cgstab:error \n";
      exit(1);
    }

    void init(double *res_in);

    scaler_cgstab(int nin, double *res_in) : base_accel() {
      n = nin;
      init(res_in);
    }

    scaler_cgstab(int nin, double *res_in, scaler_ILU &ilu) : base_accel() {
      n = nin;
      init(res_in);
      //   
      // use (LU)**(-1) res0 as res^{hat}
      //
      ilu.solve(res0, res_in);
    }

    //  destructor
    virtual ~scaler_cgstab(void);

    virtual void acc_scaler(int ia[], int ja[], double a[], int iter, 
                            double sol[], double res[], double ctol, 
                            int *iconv, const double toler[], double *rms, 
                            scaler_ILU &ilu);
                            // on entry, x^0 s.t. res^0 = b - A x^0
                            // on exit, sol = A^{-1} res^0 + x^0
}; // end scaler_cgstab


class scaler_conj: public base_accel
{
  private:
    int n;
    double  *q, *aq, *avk, *v;
    double rvkm, rmsi;

  public:

    // constructors
    scaler_conj(void) {
      cout<<"default constructor for conj grad:error \n";
      exit(1);
    }
 
    scaler_conj(int nin, double *res_in); 

    //destructor
    ~scaler_conj(void);

    virtual void acc_scaler(int ia[], int ja[], double a[], int iter, 
                            double x[], double res[], double ctol, 
                            int *iconv, const double toler[], double *rms, 
                            scaler_ILU &ilu);
}; // end  scaler_conj


// non member functions

double dot(double x1[], double x2[], int n);

void matmult(double xin[], double xout[], int n, double a[], int ia[], 
             int ja[]);

void X1EqualsX2(double* x1, double* x2, int n);

void X1EqualsX1PlusAlphaTimesX2(double* x1, double alpha, double* x2, int n);


} //end namespace SparseItObj{

#endif
