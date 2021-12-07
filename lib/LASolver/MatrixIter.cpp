#include "def_compiler.h"
#include "Standard.h"
#include "SparseItUtil.h"
#include "ILU_class.h"
#include "accel_class.h"
#include "MatrixIter.h"

using namespace SparseItObj;
using namespace std;


namespace SparseItObj { // start namespace


ItemIter::ItemIter(int val_in , ItemIter* item_ptr) 
// constructor for adding to list
{ 
  val = val_in;
  next = item_ptr;
}

ListIter::ListIter(void) // default constructor
{
  item_ptr = NULL;
  iterator = item_ptr;
}

ListIter::~ListIter(void) //destructor
{  
  ItemIter *temp_ptr, *ptr;
  ptr = item_ptr;

  while( ptr != NULL){
    temp_ptr = ptr;
    ptr = ptr->next;
    delete temp_ptr;
  }

  item_ptr = NULL;
  iterator = NULL;
}

void ListIter::iterator_start(void)
{ iterator = item_ptr; }

int ListIter::getval(void)
{
  int value = -1; // signal -tive if end of list

  if (iterator != NULL) {
    value = iterator->val;
    iterator = iterator->next;
  }

  return (value);
}

void ListIter::insert(int val)
{
  int found = 0;
  ItemIter* ptr ;
  for (ptr=item_ptr; ptr != NULL; ptr = ptr->next) {
    if (ptr->val == val) {
      found = 1;
      break;
    }
  }

  // no entry in list found, add new list entry
  if (found==0) {
    ptr = new ItemIter(val, item_ptr);
    assert(ptr != NULL);
    item_ptr = ptr;
    iterator = item_ptr;
  }
}

void ListIter::print(void)
{
  ItemIter* ptr ;
  cout << " list entries \n";
  for (ptr=item_ptr; ptr != NULL; ptr = ptr->next) {
    cout <<" list value " << ptr->val << endl;
  }
}


MatrixStruc::MatrixStruc(const int n_in, const int no_diag)
{
  ia = NULL;
  ja = NULL;
  list = NULL;

  try { // try block

    n = n_in;
    list = new ListIter [n];
    assert( list != NULL);

    stage_flag = UNPACKED;
    sort_row = 0; // flag unsorted

    if (no_diag == 0) { // insert diag entry
      int i;
      for (i=0; i<n; i++) {
        list[i].insert_fast(i); // diagonal entry in each column
        // list[i].insert(i); // diagonal entry in each column
      }
    }

  }
  catch(...) {// catch blk
    delete [] list; list = NULL;
    throw;
  }
}

MatrixStruc::~MatrixStruc(void)
{
  delete [] ia;
  delete [] ja;
  delete [] list;
}

void MatrixStruc::set_entry(int row, // global row number
                            int col // global col number
                           )
{
  if (stage_flag != UNPACKED) {
    throw General_Exception("error: data structure already compressed\n") ;
  }
  if (row < 0 || row > n-1) {
    throw  General_Exception("invalid row entry in set_entry\n") ;
  }

  if (col < 0 || col > n-1) {
    throw  General_Exception("invalid column entry in set_entry\n");
  }

  // list[col].insert(row);
  list[col].insert_fast(row); // do not chek for duplicates
}

int* MatrixStruc::getia(void) // returns new copy of ia
{
  if (stage_flag == UNPACKED) {
    pack();
  }
  int* ia_new = new int [n+1];
  assert( ia != NULL);
  int i;
  for (i=0; i<=n; i++) {ia_new[i] = ia[i];}

  return ( ia_new );
}

int* MatrixStruc::getja(void) // returns new copy of ja
{
  if (stage_flag == UNPACKED) {
    pack();
  }
  int nja = ia[n];
  int* ja_new = new int [nja];
  assert( ja_new != NULL);
  int i;
  for (i=0; i<nja; i++) {ja_new[i] = ja[i];}

  return (ja_new);
}

int MatrixStruc::getnja(void)
{
  if (stage_flag == UNPACKED) {
    pack();
  }
  return (ia[n]);
}

void MatrixStruc::pack(void)
{
  if (stage_flag == PACKED) {
    throw General_Exception("error: data structure already packed\n");
  }

  int* count  = NULL;
  count = new int [n];
  assert( count != NULL);

  ia = NULL;
  ja = NULL;

  try { // try blk

    int i;
    for (i=0; i<n; i++) {
      count[i] = 0;      // count of nonzeros in each row
    }// end loop

    for (i=0; i<n; i++) { // scan over columns
      list[i].iterator_start(); 
      int irow;
      irow = list[i].getval();
      while (irow != -1) {
        count[irow]++;
        irow = list[i].getval();
      }
    }

    ia = new int [n+1];
    assert(ia != NULL);

    ia[0] = 0;
    for (i=0; i<n; i++) {
      ia[i+1] = ia[i] + count[i];
    }

    ja = new int [ia[n]];
    assert(ja != NULL);

    for (i=0; i<n; i++) {count[i] = ia[i+1]-1;}

    for (i=n-1; i>=0; i--) {// scan over list of columns in reverse order
      list[i].iterator_start();
      int irow = list[i].getval();
      while (irow != -1) {
        assert(count[irow] >= ia[irow]);

        ja[count[irow]] = i;
        count[irow]--;
        irow = list[i].getval();
      }
    }

    // ia,ja list may contain duplicates
    // but  we can remove these
    MatrixIter_ia_ja_remove_duplicates(ia, ja, n);

    // rows now sorted
    sort_row = 1;

  } // end try blk
  catch(...) { // catch anything
    delete [] count;
    count = NULL;
    delete [] ia; ia = NULL;
    delete [] ja; ja = NULL;
    delete [] list; list = NULL;
    throw;
  } // end catch

  delete [] list;
  list = NULL;

  delete [] count;

  stage_flag = PACKED;
}

double& MatrixIter::aValue(const int row, // global row number of entry
                           const int col  // global col number of entry
                          )
{
  int i;
  int found = -1;
  for (i=ia[row]; i<ia[row+1]; i++) { // start loop
    if (ja[i] == col) {
      found = i;
      break;
    }
  }// end loop
  if (found == -1) {
    throw General_Exception("error: attempt to insert entry not in structure\n");
  }

  return( a[ found ] );
}

double& MatrixIter::aValue_bsearch(const int row, // global row number
                                   const int col  // global col number
                                  ) // uses binary search to insert
                      // probably only necessary if large number of entries
                      // (>20) in a row
{
  if (sort_row == 0) {
    throw General_Exception("error: binary search cannot be used rows unsorted\n");
  }
          
  int size = ia[row+1] - ia[row];

  int found = bsearchIter(&(ja[ia[row]]), size, col);

  if (found == -1) {
    throw General_Exception("error: attempt to insert entry not in structure\n");
  }

  return (a[found]);
}

int MatrixIter::check_entry(int row, // row index
                            int col  // col index
                           )
{
  int i;
  int found = -1;
  for (i=ia[row]; i<ia[row+1]; i++) { // start loop
    if (ja[i] == col){
      found = i;
      break;
    }
  }// end loop

  if (found == -1) { 
    return 0;
  }
  else{
    return +1;
  }
} 

MatrixIter::MatrixIter(MatrixStruc& iaja_set) // iaja data structure
{
  sort_row = iaja_set.get_sort_row();

  int n_in = iaja_set.getn();
  int* ia_temp = NULL;
  int* ja_temp = NULL;

  try{
    ia_temp = iaja_set.getia();
    ja_temp = iaja_set.getja();
    init(n_in, ia_temp, ja_temp);
  }
  catch(...) {
    delete [] ia_temp; ia_temp = NULL;
    delete [] ja_temp; ja_temp = NULL;
    throw;
  }

  delete [] ia_temp;
  delete [] ja_temp; // clean up temp ia, ja
                     // copies in MatrixIter class now
}

MatrixIter::MatrixIter(const int n_in, const int* ia_in, const int* ja_in)
{
  sort_row = 0; // flag as possibly unsorted
  init(n_in, ia_in, ja_in);
}

void MatrixIter::init( const int n_in, const int* ia_in, const int* ja_in)
{
  a = NULL;
  b = NULL;
  res = NULL;
  toler = NULL;
  a_scal = NULL;
  ia = NULL;
  ja = NULL;
  ilu_ptr = NULL;
  lorder_user = NULL;
  porder_user = NULL;

  try { // try blk
      
    n = n_in;
    ia = new int [n+1];
    assert(ia != NULL);

    ja = new int [ia_in[n]];
    assert(ja != NULL);

    int i;
    for (i=0; i<= n; i++) {
      ia[i] = ia_in[i];
    }
    for (i=0; i<= ia_in[n]-1; i++) {
      ja[i] = ja_in[i];
    }

    sym_factor_status = UNFACTORED;
    num_factor_status = UNFACTORED;

    a = new double [ia[n]];
    assert(a != NULL);
    b = new double [n];
    assert(b != NULL);

    zeroa();
    zerob();

    a_scal = new double [n];
    assert(a_scal != NULL);

    toler = new double [n];
    assert(toler != NULL);

    { int i;
      for (i=0; i<=n-1; i++) {
        toler[i] = 0.0;
        a_scal[i] = 1.0;
      }
    }

    res = new double [n];
    assert( res != NULL);

    nzero = 0;

    lorder_user = NULL;  // flag use of internal ordering
    porder_user = NULL;

  }// end try blk
  catch(...) { // catch
    delete [] a; a = NULL;
    delete [] b; b = NULL;
    delete [] toler; toler = NULL;
    delete [] res; res = NULL;
    delete [] a_scal; a_scal = NULL;
    delete [] ia; ia = NULL;
    delete [] ja; ja = NULL;
    delete ilu_ptr; ilu_ptr = NULL;
    delete [] lorder_user; lorder_user = NULL;
    delete [] porder_user; porder_user = NULL;
    throw;
  }// end catch

}// end init()

MatrixIter::~MatrixIter(void)
{
  delete [] a; a = NULL;
  delete [] b; b = NULL;
  delete [] toler; toler = NULL;
  delete [] res; res = NULL;
  delete [] a_scal; a_scal = NULL;
  delete [] ia; ia = NULL;
  delete [] ja; ja = NULL;
  delete ilu_ptr; ilu_ptr = NULL;
  delete [] lorder_user; lorder_user = NULL;
  delete [] porder_user; porder_user = NULL;
}

void MatrixIter::set_toler(const double* tol_in) // tolerance vector
{
  // set convergence tolerance for inner
  // iteration relative to size of variable

  int i;

  for (i=0; i<=n-1; i++) {
    toler[i] = tol_in[i];
  }
}

void MatrixIter::sfac(ParamIter& param) //, ofstream& oFile)
{
  if (ilu_ptr != NULL) {
    delete ilu_ptr;
    ilu_ptr = NULL;
  }
      
  { int* lorder = NULL; // start ordering block
    lorder = get_lorder( param );
    try {
      ilu_ptr = new scaler_ILU(n, lorder, porder_user);
    }
    catch(...) {
      delete [] lorder;
      lorder = NULL;
      throw;
    }
    assert(ilu_ptr != NULL);
    delete [] lorder;
  } // end ordering block

  int ier = 0;

  ilu_ptr->sfac2(ia, ja, param.level, &nzero, &ier);
  if (ier != 0) {
    throw General_Exception("error in symbolic level ILU \n");
  }

  // if (param.info != 0) {
  //   oFile << " nonzeros in factors: " << nzero << endl;
  // }

  sym_factor_status = FACTORED;
  num_factor_status = UNFACTORED;
}

int* MatrixIter::get_lorder(ParamIter& param)
{
  int* lorder;

  lorder = new int [n];
  assert (lorder != NULL);

  if (lorder_user != NULL) { // use user ordering
    int i;
    for (i=0; i<=n-1; i++) {
      lorder[i] = lorder_user[i];
    }

  }
  else{ // use internal ordering
    if (param.order == 0){ // natural ordering
      int i;
      for (i=0; i<=n-1; i++) {
        lorder[i] = i;
      }
    }
    else{ // rcm odering
      rcm(ia, ja, n, lorder);
    }
  }

  return lorder;
}

void MatrixIter::set_user_ordering_vector_column(int* porder_in)
             // user ordering vector porder[new_order] = old_order
             // lorder - n-length array allocated and deallocated
             // by user. Overrides any other ordering option
             // if set  /
             // column ordering (normall col = row by default)
{
  if (porder_in == NULL) {
    return;
  }

  porder_user = new int [n];
  assert(porder_user != NULL);

  int i;
  for (i=0; i<n; i++) {
    porder_user[i] = porder_in[i];
  }

  // chek valid ordering
  { int* ichek = NULL;
    ichek = new int [n];
    assert( ichek != NULL);
    int ii;
    for (ii=0; ii<n; ii++) { // start loop over nodes
      ichek[ii] = -1;
    }// end loop over nodes
    for (ii=0; ii<n; ii++) { // start loop over nodes
      int jj = porder_user[ii];

      if (jj < 0 || jj >= n) {
        throw  General_Exception("error_1 in user  ordering\n");
      }

      if (ichek[jj] != -1) {
        throw  General_Exception("error_2 in user  ordering\n");
      }
      ichek[jj] = +1;

    }// end loop over nodes

    delete [] ichek;
  }

}

void MatrixIter::set_user_ordering_vector(int* lorder_in)
             // user ordering vector lorder[new_order] = old_order
             // lorder - n-length array allocated and deallocated
             // by user. Overrides any other ordering option
             // if set

{
  if (lorder_in == NULL) {
    return;
  }

  lorder_user = new int [n];
  assert(lorder_user != NULL);

  int i;
  for (i=0; i<n; i++) {
    lorder_user[i] = lorder_in[i];
  }

  // chek valid ordering
  { int* ichek;
    ichek = new int [n];
    assert(ichek != NULL);
    int ii;
    for (ii=0; ii<n; ii++) { // start loop over nodes
      ichek[ii] = -1;
    }
    for (ii=0; ii<n; ii++) { // start loop over nodes
      int jj = lorder_user[ii];

      if (jj < 0 || jj >= n) {
        throw  General_Exception("error_3 in user  ordering\n");
      }

      if (ichek[jj] != -1) {
        throw  General_Exception("error_4 in user  ordering\n");
      }
      ichek[jj] = +1;

    }// end loop over nodes

    delete [] ichek;
  }

}

void MatrixIter::solveWithOldFactors(ParamIter& param, // solution parameters
                                                       // allocated in caller
                                     double* x, // solution (n-length array)
                                                // allocated in caller
                                     int &nitr, // return number of iterations
                                                // return -1 if not converged
                                    //  ofstream& oFile, 
                                     const int initial_guess
                  // = 0 assume x = init guess = 0
                  //    NOTE: x is zeroed in solve for this option
                  //    returns x = A^{-1}b
                  // !=0 assume x is intial guess
                  //    returns x = A^{-1}res^0 + x^0
                  //    res^0 = b - A x^0
                                   )
{
  num_factor_status = FACTORED; // use old factors

  // solve(param, x, nitr, oFile, initial_guess); 
  solve(param, x, nitr, initial_guess); 
     
}

void MatrixIter::solve(ParamIter& param,
               double* x, // solution
               int &nitr, // return number of iterations
                          // return -1 if not converged
              //  ofstream& oFile, 
               const int initial_guess
                  // = 0 assume x = init guess = 0
                  //    NOTE: x is zeroed in solve for this option
                  //    returns x = A^{-1}b
                  // !=0 assume x is intial guess
                  //    returns x = A^{-1}res^0 + x^0 
                  //    res^0 = b - A x^0
                      )
{
  if (num_factor_status == FACTORED) { // solve with old factors
    int i;
    for (i=0; i<n; i++) {
      b[i] *= a_scal[i];
    }
    goto OLD_FACTORS;
  }

  if (param.drop_ilu == 0 && sym_factor_status != FACTORED) {
    throw  General_Exception("error: solve called with no symbolic ILU \n");
  }

  if (param.iaccel == -1) {
    param.iscal = 0;  // no scaling for conjugate grad
                      // could use symm scaling, but not
                      // bothering now
  }

  if (param.iscal != 0) {// start scaling if
    if (param.ipiv == 0) {// piv if
      scal( n, ia, ja, a, b, a_scal,0);
    }
    else {
      scal(n, ia, ja, a, b, a_scal, 1);
    }// end piv if
  }// end scal
  else {
    int i;
    for (i=0; i<=n-1; i++) {
      a_scal[i] = 1.0;
    }
  }

  if (param.drop_ilu != 0) { // begin drop tol stuff
    if (sym_factor_status == UNFACTORED) {// redo drop tol ilu
      if (ilu_ptr != NULL) {delete ilu_ptr; ilu_ptr = NULL;}

      { int* lorder; // start ordering block
        lorder = get_lorder(param);
        ilu_ptr = new scaler_ILU(n, lorder, porder_user);
        assert(ilu_ptr != NULL);
        delete [] lorder;
      } // end ordering block

      ilu_ptr->setDropType(3); // drop relative to current row

      if (param.ipiv == 1) {
        ilu_ptr->setpivoton();
      }

      int ier = 0;
      ilu_ptr->facdrp2(ia, ja, a, param.drop_tol, &nzero, &ier);

      if (ier != 0) {
        throw General_Exception("error in drop tol ILU\n");
      }

      sym_factor_status = FACTORED;
      num_factor_status = FACTORED;

    }// end redo drop tol ilu
  }// end drop tol stuff

  if (num_factor_status == UNFACTORED) {
    ilu_ptr->factor(ia, ja, a);  // numeric factor
    num_factor_status = FACTORED;
  }

  // iterations

  OLD_FACTORS: ;  // skip to here if using old factors

  { // acceleration block

    // zero intial guess
    double rmsi = 0.0;

    if (initial_guess == 0) { // x^0 = 0

      { int i;// start block
        for (i=0; i<=n-1; i++) { // start loop
          x[i] = 0.0;
                 // for nonzero initial guess
                 // compute res = b - A x
          res[i] = b[i];
          rmsi += res[i]*res[i];
        }
      }
    } 
    else { // x^0 is initial guess != 0
      matmult(x, res, n, a, ia, ja); // res = Ax
      int i;
      for (i=0; i<n; i++) {// start loop
        res[i] = b[i] - res[i]; // form res^0 = b - Ax^0
        rmsi += res[i]*res[i];
      }
    }

    rmsi = sqrt( rmsi);

    // scaler_cgstab* accel_ptr;
    base_accel* accel_ptr;

    switch(param.iaccel) { // start switch

      case 0: { // cgstab acceleration
        if (param.new_rhat == 0) {
          accel_ptr = new scaler_cgstab(n, res);
        }
        else {
          accel_ptr = new scaler_cgstab(n, res, *ilu_ptr); 
          //new scaler_cgstab(n, res, *ilu_ptr);
          // uses(LU)^{-1}*b as rhat in cgstab acceleration
          // avoids possible problems with initial guess
        }
        break;
      }

      case 1: { // orthomin 
        accel_ptr = new scaler_orthomin(n, res, param.north);
        break;
      }

      case -1: { // conj gradient
        accel_ptr = new scaler_conj(n, res);
        break;
      }

      default: {
        throw  General_Exception("error: invalid param.iaccel\n");
      }

    }// end switch

    // oFile.setf(ios::left, ios::adjustfield); // left justify
    // oFile.precision(5);

    // if (param.info == 1) {
    //   oFile << "initial   rms: "  << rmsi << endl;
    // }

    nitr = 0;

    { // iteration block
      double rms;
      int iter;
      int iconv;

      for (iter=1; iter<= param.nitmax; iter++) {// iteration loop
        nitr++;

        accel_ptr->acc_scaler(ia, ja, a, iter, x, res, param.resid_reduc, 
                              &iconv, toler, &rms, *ilu_ptr);

        // if (param.info == 1) {
        //   oFile << "iter: " << iter << "  " << " rms: " << rms << endl;
        // }

        if (iconv == 0) {break;}
      }

      if (iconv != 0) {nitr = -1;} // flag non-convergence

    }// end iteration block

    delete accel_ptr;
  } 

  num_factor_status = UNFACTORED; // force new numeric factor
                                  // at next call
}

void MatrixIter::set_row(const int i, // set values of i'th row
                         double* row  // row is full storage mode
                                      // for i'th row
                        )
{ 
  int j;
  for (j=ia[i]; j<= ia[i+1]-1; j++) {
    int id = ja[j];
    a[j] = row[ id ];
  }
}
 
void MatrixIter::zero_row(const int i, // i'th row
                          double* row  // full storage n-length
                                       // array.  Zero those entries
                                       // in row corresponding to nonzeros
                                       // in row i of matrix
                         )
{
  int j;
  for (j=ia[i]; j<= ia[i+1]-1; j++) {
    int id = ja[j];
    row[ id ] = 0.0;
  }
}

int MatrixIter::get_upper_fill(void)
{return ilu_ptr->UpperFill;}

int MatrixIter::get_lower_fill(void) 
{return ilu_ptr->LowerFill;}

void MatrixIter::zeroa(void)
{
  int nja = ia[n];
  int i;
  for (i=0; i<nja; i++) {a[i] = 0.0;}
}

void MatrixIter::zerob(void)
{ 
  int i;
  for (i=0; i<n; i++) {b[i] = 0.0;}
}

double MatrixIter::mult_row(const int row, // input row number
                            double* val    // array of size [n] with values
                           )
{
  if ((row < 0) || (row >= n)) {
    throw General_Exception("error in MatrixDirect::mult_row\n");
  }

  double result = 0.;

  for (int j=ia[row]; j<= ia[row+1]-1; j++) {
    int id = ja[j];
    result += a[j]*val[id];
  }

  return( result );
}

void ListIter::insert_fast(int val)
{
  //insert into list, do not chek for duplicates
  ItemIter* ptr;

  ptr = new ItemIter(val, item_ptr);
  assert(ptr != NULL);
  item_ptr = ptr;
  iterator = item_ptr;
}

void MatrixIter_ia_ja_remove_duplicates(int* & ia, int* & ja, int n)
//
//  on exit, ia,ja replaced by sorted lists
// with duplicates removed
{
  int* count = new int[n];
  assert(count != NULL);

  int* ia_new = NULL;
  int* ja_new = NULL;

  try { // start block

    {// start block

      int i;

      for (i=0; i<n; i++) {count[i] = 0;}

      for (i=0; i<n; i++) { // loop over nodes

        // list should be sorted, but chek 
        { // check block
          int k;
          for (k=ia[i]; k<ia[i+1]; k++) { //loop over neighs
            if (k > ia[i]) {
              assert( ja[k] >= ja[k-1]);
            }
          }
        }

        // shell( &ja[ ia[i] ] , (ia[i+1] - ia[i] ) );
        int j;
        for (j=ia[i]; j<ia[i+1]; j++) { //loop over neighs
          int test_ok = 1;

          if (j+1 < ia[i+1]) {
            if (ja[j] == ja[j+1]) {
              test_ok = 0;
            }
          }

          if (test_ok == 1) {
            count[i]++;
          } 
        }
      }

      ia_new = new int [n+1];
      assert(ia_new != NULL);

      ia_new[0] = 0;
      for (i=0; i<n; i++) {
        ia_new[i+1] = ia_new[i] + count[i];
      }

      ja_new = new int [ia_new[n]];
      assert(ja_new != NULL);

      for (i=0; i<n; i++) {count[i] = 0;}

      for (i=0; i<n; i++) { // loop over nodes
        int j;
        for (j=ia[i]; j<ia[i+1]; j++) { //loop over neighs
          int test_ok = 1;

          if (j+1 < ia[i+1]) {
            if (ja[j] == ja[j+1]) {
              test_ok = 0;
            }
          }

          if (test_ok == 1) {
            assert(ia_new[i]+count[i] < ia_new[i+1] );

            ja_new[ ia_new[i]+count[i] ] = ja[j];
            count[i]++;
          }
        }
      }
    }
  }
  catch(...) {// start catch
    delete [] ia_new; ia_new = NULL;
    delete [] ja_new; ja_new = NULL;
    delete [] count; count = NULL;
    throw;
  }

  delete [] ia;
  delete [] ja;

  ia = ia_new;
  ja = ja_new;

  delete [] count;
}

} // end namespace SparseItObj{
