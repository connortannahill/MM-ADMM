
#ifndef MATRIX_ITER_INC
#define MATRIX_ITER_INC

#include <iostream>

namespace SparseItObj{ // start namespace SparseitObj

void rcm(const int ia[], const int ja[], int n, int lorder []);


class ItemIter {

  friend class ListIter;

  private:
    int val;
    ItemIter *next;

    // default constructor
    ItemIter(void) {
      std::cout <<"error: default constructor for ItemIter called\n";
      exit(1);
    }

    // constructor for adding to list
    ItemIter(int val_in, ItemIter* item_ptr); 

    // copy constructor
    ItemIter(const ItemIter &item_in) {
      std::cout << "error: copy constructor for ItemIter called\n";
      exit(1);
    }

    // destructor
    ~ItemIter(void){ }

};


class ListIter{
 
  private:
    ItemIter *item_ptr;

    ItemIter *iterator;

  public:
    ListIter(void); // default constructor
    ~ListIter(void);  //destructor
 
    void iterator_start(void);

    int getval(void);

    // insert with no duplicates
    void insert(int val);

    // insert, do not chek for duplicates
    void insert_fast(int val);

    void print(void);
  };


class MatrixStruc{

  private:

    int  n; // number of nodes
    int* ja;
    int* ia;
    ListIter* list;  // linked list of
                     // row entries in each column

    int  sort_row; // = 0 not sorted
                   // = 1 sorted

    enum Stage {UNPACKED, PACKED};
    Stage stage_flag;

    MatrixStruc(const MatrixStruc &); // no copy constuctor
    MatrixStruc(void); // no def constructor

  public:

    MatrixStruc(const int n_in, // number of unknowns
                                // i.e. A is an n x n matrix
                const int no_diag = 0  // = 0 add nozeros to diag
                                       // = 1 do not
    );

    ~MatrixStruc(void);

    int get_sort_row(void) {return sort_row;}

    void set_entry(int row, // global row number, 0<= row <= n-1
                   int col  // global col number, 0<= col <= n-1
                  ); //set (i,j) as nonzero entry in matrix

    int* getia(void); // returns new copy of ia 
    int* getja(void); // returns new copy of ja
    int getnja(void); // return nonzeros in a() - 1
    void pack(void);  // convert from linked list -> compressed row
    int getn(void) const { return n;}

}; // end MatrixStruc


class scaler_ILU; // forward declaration


class ParamIter{

  public:

    int order; // = 0 natural
               // = 1 RCM

    int level; // level of ilu

    int drop_ilu; // = 0 level of fill ilu
                  // = 1 drop tol
   
    int ipiv; // = 0 not pivoting
              //  = 1 pivoting if drop_ilu=1

    int iscal; // = 0 no scaling 
               // = 1 scaling

    int nitmax; // max number of inner itns

    double resid_reduc; // residual reduction convergence criteria

    int info; // = 0 no inner itn info
              // = 1 inner itn info written to "output"

    double drop_tol; // drop tolerance 

    int new_rhat; // = 0 use rhat = r^0
                  // = 1 use rhat = (LU)^{-1} r^0
                  // vector in cgstab algorithm
                  // i.e. cgstab residuals orthogonal
                  // to Krylov space P(A^t) rhat

    int iaccel; // = 0 cgstab (default)
                // = 1 orthomin
                // = -1 conj gradient

    int north; // number of orthogs for
               // orthomin

    // default constructor
    ParamIter( void) {
      order = 1;
      level = 1;
      drop_ilu = 0;
      iscal = 1;
      nitmax = 30;
      resid_reduc = 1.e-6;
      drop_tol = 1.e-3;
      info = 1;
      new_rhat = 0;
      iaccel = 0;
      north = 10;
      ipiv = 0;
    } // end default constructor

    private:
      ParamIter(const ParamIter &); // no copy constuctor

    public:
      ~ParamIter(void) { } // destructor

}; // end class param


class MatrixIter {

  private:
    int n; // number of unknowns
    enum stageType{FACTORED, UNFACTORED};
    stageType sym_factor_status; // symbolic factor status
    stageType num_factor_status; // numeric factor status

    int sort_row; // = 0 not sorted
                  // = 1 sorted

    double* a; // matrix
    double* b; // rhs

    double* res; // residual of linear eqns
    double* toler;  // tolerance vector
    double* a_scal; // scaling vector

    int* ia; 
    int* ja; // ia,ja arrays for the matrix
             // in gen'l can be different from
             // grid ia-ja

    scaler_ILU* ilu_ptr; // ptr to ILU type
                         // defined in "ILU_class.h"

    int nzero; // nonzeros in the ILU

    int* lorder_user; // user input ordering vector
                      // lorder[new_order] = old_order
                      // this is actually the
                      // row ordering, by default
                      // the same ordering is used
                      // for both rows and columns
                      // (symmetric reordering0
    int* porder_user; // col ordering
                      // used only if col pivoting used
  public:

    MatrixIter(const MatrixIter &) {
      std::cout << "error copy constructor  MatrixIter  called\n";
      exit(1);
    }

    //  default constructor
    MatrixIter(void) {
      std::cout << "error default constructor  MatrixIter called\n";
      exit(1);
    }

    MatrixIter(MatrixStruc& iaja_set); // iaja data structure

    void init(const int n_in, const int* ia_in, const int* ja_in);

    MatrixIter(const int n_in, const int* ia_in, const int* ja_in);

    // destructor
    ~MatrixIter(void);

    // ilu functions
    void sfac(ParamIter& param//, // ILU parameters
              //ofstream& oFile
              ); // symbolic factor for level of fill ILU
  
    void solve(ParamIter& param, // solution parameters
                                 // allocated in caller
               double* x, // solution (n-length array
                          // allocated in caller
               int &nitr, // return number of iterations
                          // return -1 if not converged
               //ofstream& oFile, 
               const int initial_guess = 0
                  // = 0 assume x = init guess = 0
                  //     NOTE: x is zeroed in solve for this option
                  //     returns x = A^{-1}b
                  // !=0 assume x is intial guess
                  //     returns x = A^{-1}res^0 + x^0
                  //     res^0 = b - A x^0
              );
    
    void printRHS(int n);
    void printMat(int n);

    void solveWithOldFactors(ParamIter& param, // solution parameters
                                               // allocated in caller
               double* x, // solution (n-length array
                          // allocated in caller
               int &nitr, // return number of iterations
                          // return -1 if not converged
              //  ofstream& oFile,
               const int initial_guess = 0
                  // = 0 assume x = init guess = 0
                  //     NOTE: x is zeroed in solve for this option
                  //     returns x = A^{-1}b
                  // !=0 assume x is intial guess
                  //     returns x = A^{-1}res^0 + x^0
                  //     res^0 = b - A x^0
              );

    void set_user_ordering_vector(int* lorder);
                 // user ordering vector lorder[new_order] = old_order
                 // lorder - n-length array allocated and deallocated
                 // by user. Overrides any other ordering option
                 // if set

    void set_user_ordering_vector_column( int* porder_in);
                 // user ordering vector porder[new_order] = old_order
                 // lorder - n-length array allocated and deallocated
                 // by user. Overrides any other ordering option
                 // if set 
                 // column ordering (normall col = row by default)

    void set_toler(const double* tol_in); // tolerance vector
                                          // n-length array allocated in caller

    void set_row(const int i, // set values of i'th row
                 double* row  // row is full storage mode
                              // for i'th row (i.e. an n=length 
                              // vector)
                );

    void zero_row(const int i, // i'th row
                  double* row  // full storage n-length
                               // array.  Zero those entries
                               // in row corresponding to nonzeros
                               // in row i of matrix
                 );

    int* get_lorder(ParamIter& param);

    double& bValue(const int i)
    {return b[i];}  // access to right hand side

    double& aValue(const int k)
    {return a[k];}  // int k determined with the help
                    // of rowBegin, rowEndPlusOne,
                    // getColIndex  functions
                    // efficient access to entries in sparse
                    // matrix
 
    double& aValue(const int row, // global row number of entry
                   const int col  // global col number of entry
                  ); // access to entry in sparse matrix

    double& aValue_bsearch(const int row, // global row number
                           const int col  // global col number
                          ); // uses binary search to insert
                             // probably only necessary if large number 
                             // of entries (>20) in a row
  
    int check_entry(int i , // row index
                    int j // col index
                   ); // returns 0 if requested (i,j) pair not in data stucture
                      // returns !0 if ok

    void zeroa(void); // set all entries in A to zero

    void zerob(void); // seta all entries in rhs to zero

    void set_sym_unfactored(void)
    {
      sym_factor_status = UNFACTORED;

      // for drop tol ilu, force new drop tol ilu,
      // otherwise, uses previous sparsity pattern
    }  

    int get_sym_factor_status(void) {
      if (sym_factor_status == FACTORED) {return (1);}
      else {return (0);}
    }

    int getnonzero (void) {return nzero;}

    int get_upper_fill(void);

    int get_lower_fill(void); 

    // the following functions are used to access and set
    // sparse matrix entries;
    // In row i, entries in are stored in a[k], such
    // that rowBegin(i) < = k < rowEndPlusOne(i);
    // So, for MatrixIter matrix, if a_{i,j} corresponds
    // to a nonzero in row i and column j, then:
    //
    //       matrix.aValue(k) = a_{i,j}, with j = matrix.getColIndex(k)
    //
    //       matrix.rowBegin(i) < = k < matrix.rowEndPlusOne( i )

    int rowBegin(const int row) const {return ia[row];}

    int rowEndPlusOne(const int row) const {return ia[row+1];}
   
    int getColIndex (const int k) const {return ja[k];}

    double mult_row(const int row, // input row number
                    double* val    // array of size [n] with values
                   );

    int get_n(void) {return (n);}

    int* get_ia(void) {return (ia);}

    int* get_ja(void) {return (ja);}

};// end class MatrixIter


void MatrixIter_ia_ja_remove_duplicates(int* & ia, int* & ja, int n);

}// end namespace SparseItObj

#endif
