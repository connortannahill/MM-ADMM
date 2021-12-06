
#ifndef ILU_CLASS_INC
#define ILU_CLASS_INC

// change to float for single precision
// factors for scaler ilu

namespace SparseItObj { // start namespace SparseItObj

typedef double SizeScalerFactor;

/*
    each row is held as a structure
    so, for i'th row, with "nz" nonzero entries
 
       for k=0,.., nz-1
       jaf[k] = col index
       lev[k] = level of entry
       af[k]  = value
       k = diag is the location of the diagonal
        of the row  ie  af[ row.diag ] = diag entry
*/

class base_sparse_row
{
  public:

    int *jaf;
    int *lev;
    int nz;
    int diag;

    // default constructor
    base_sparse_row(void) {
      jaf = NULL;
      lev = NULL;
      nz = 0;
      diag = -1;
    }

    //   default destructor
    virtual ~base_sparse_row(void) {
      if (jaf != NULL) {
        delete [] jaf;
        jaf = NULL;
      }
      if (lev != NULL) {
        delete [] lev;
        lev = NULL;
      }
      nz = 0;
      diag = -1;
    }

  private:

  // default copy constructor should not be called
  base_sparse_row(const base_sparse_row &row);

};  // end base_sparse_row


class scaler_sparse_row:public base_sparse_row
{
  public:

    SizeScalerFactor *af;

    // default constructor
    scaler_sparse_row(void):base_sparse_row()
    {
      af = NULL;
    }
 

  private:
  // default copy constructor should not be called
  scaler_sparse_row(const scaler_sparse_row &row);

  public:
 
  // destructor
    ~scaler_sparse_row(void) {
      if (af != NULL) {
        delete [] af;
        af = NULL;
      }
    }

}; // end scaler_sparse_row


class base_ILU
{
  protected:

    int *lorder, *invord; // row ordering

    int* porder;  // porder[new]  old   col ordering
    int* invpord; // invpord[old] = new
    int n;
 
    enum space {allocated, unallocated };
    space Dspace;
 
    enum symfac {symfac_done, symfac_not_done};
    symfac Symfac;
 
    enum numfac {numfac_done, numfac_not_done};
    numfac Numfac;

  public:

    int LowerFill;
    int UpperFill;

  private:

    base_ILU(void); // def constructor should not be called
 
  public:
 
    base_ILU(int n_in, int lord_in[], int* pord_in);
 
    // destructor
    ~base_ILU(void){ 
      delete [] lorder;
      delete [] invord; 
      delete [] porder;
      delete [] invpord;
    }

    int check_factor(void) {
      if (Numfac == numfac_not_done) {
        exit(1);
      }
      else {
        return 0;
      }
    }
 
}; // end class base_ILU
 

class scaler_ILU:public base_ILU
{
  protected: 

  // ptr to structure of ILU
  scaler_sparse_row *rowsp;

  int DropType; // = 0 sqrt  (diag (orig row) diag(orig col) )
                // = 1 max orig row
                // = 2 orig diag
                // = 3 current diag
                  
  int PivType; // =0 no piviting
               // =1 pivoting

  public:

    // constructors

  private:

    scaler_ILU(void); // should not be called

  public:

    void init(void);

    scaler_ILU(int n_in, int lord_in[], int* pord_in)
     : base_ILU( n_in, lord_in,  pord_in)
    {  
      rowsp = NULL;
      init();
    }

    // destructor
    void clean(void);
    ~scaler_ILU(void) {clean();}

  // member functions
  void setDropType(int drop_in) // see DropType
  {
    DropType = drop_in;
    assert(DropType >= 0 || DropType < 4);
  }

  void setpivoton(void)
  {PivType = 1;}

  void setpivotoff(void)
  {PivType = 0;}

  // merge2 called by sfac2
  void merge2(int i, int list[], int first, int lrow[], int level);

  // symbolic factorization
  void sfac2(const int ia[], const int ja[], int level, int *nzero, int *ier);

  // numeric factorization
  void factor(const int ia[], const int ja[], const double a[]);

  // solve with ILU factors
  void solve( double x[], const double b[]);  

  void merged(int i,  int list[], int first, double arow[], double dropt,
              double odiag[]); // used by facdrp
                               // obsolete

  void facdrp(const int ia[], const int ja[], double a[], double dropt, 
              int *nzero, int *ier); // this is obsolete

  void facdrp2(const int ia[], const int ja[], double a[], double dropt, int 
               *nzero, int *ier); // new drop tol factor, allows for pivoting

  void elimrow(const int i ,   // row being elimiated
               double* odiag,  // dtop tol test data
               double* arow,   // full length array of
                               // nonzeros for this row
               int* marker,    // nonzero in col_j -> marker[j] = i
               int& count_up,
               int& count_low, // col entries for this row
               int* lower,     // for col < diag,
               int* upper,     // lower[0,..,count_low]
                               // for col >=diag
                               // upper[0,..,count_up]
                               // on entry, undered
                               // on exit, lower[] sorted
                               // upper[] unsorted
               double dropt    // drop tol
              ); // used by facdrp2
 
}; // end class scaler_ILU


void scal(const int n, // number of unknowns
          int* ia, int* ja, 
          double* a,        // orig matrix
          double* b,        // right hand side
          double* scal_fac, // scaling factors
          int scal_method   // = 0 scal by diag
                            // = 1 scal by max row
         );

} //end namespace SparseItObj

#endif
