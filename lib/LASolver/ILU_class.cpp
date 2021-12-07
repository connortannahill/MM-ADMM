#include "Standard.h"
#include "SparseItUtil.h"
#include "ILU_class.h"
#include "def_compiler.h"

using namespace SparseItObj;
using namespace std;;


/*    function definitions     */


//   scaler class functions

namespace SparseItObj{ // start namespace

void scaler_ILU::merge2(int i, int list[],int first, int lrow[], int level)
/*
  input
    
    i       current row being factored (Dolittle form)
    list[]  implied linked list of nonzeros in current
            row.  On entry, only original nonzero columsn
            in this row.  1st nonzero column = fisrst
            , list[fisrt], list[ list[first]], ...., until
            list[  [....] ] = n+1 -> signals end of list
    n       number of rows
    first   (see list[]) fist nonzero in this row
    lrow[]  level of fill for entries in this row, nonzeros
            correspond to those in list[]
            Full storage vector (size = n)
            All original nonzeros have lrow = 0
            Otherwise, lrow = infinity
    level   max level of fill
    *rowsp  ptr to global structure containing info about
            previous rows, and level of fill of entries
            in prev rows

  output

    list[]   now contains all nonzeros in this factored row i
             (implied linked list)
*/
{
  int row, next, oldlst,  ii;
  int nxtlst, levup, levlow, levnew;


  next = first;

  while (next < i) {

    oldlst = next;
    nxtlst = list[next];
    row = next;
 
    levlow = lrow[ oldlst ];

    /* scan row "row" of U */
 
    for (ii=rowsp[row].diag+1; ii<=rowsp[row].nz-1; ii++) {
 
      /* skip thru linked list until we find possible insertion pt */
       while (rowsp[row].jaf[ii] > nxtlst) {
         oldlst = nxtlst;
         nxtlst = list[oldlst]; 
       }
                
       if (rowsp[row].jaf[ii] < nxtlst) {
         levnew = levlow + rowsp[row].lev[ii] + 1;
         if (levnew <= level) {
           list[ oldlst] = rowsp[row].jaf[ii];
           list[ rowsp[row].jaf[ii]  ] = nxtlst;
           oldlst = rowsp[row].jaf[ii];
           lrow[ oldlst ] = levnew;
         }
       }
       else {
         oldlst = nxtlst;
         levup = rowsp[row].lev[ii];

         /* take min of old and new levels */

         lrow[ oldlst ] = levup + levlow + 1 < lrow[oldlst] ?
                          levup + levlow + 1 : lrow[oldlst];
         nxtlst = list[oldlst ];
       }
     }
     next = list[ next ];
  }

  return  ;
}


/* symblic factor */

void scaler_ILU::sfac2(const int ia[], const int ja[], int level, 
int *nzero, int *ier)

/*   level based symbolic factor

  input

    ia, ja     structure of a
    n          number of unknowns (class variable)
    lorder[]   lorder[new_order] = old_order (class varaible)
    invord[]   invord[old_order] = new_order (class variable)
                    row ordering
    porder[]   porder[ new_order] = old_order
    invpord[]  invpord[old_order] = new_order
                    col ordering

    level      max level of fill
  
  output:

    *nzero     number of nonzeros in ILU
    *ier       error flag
    *rowsp     stucture contains symbolic ILU info
               (class variable)
    Symfac     (class variable) indicates symfac done ok
    Dspace     (class variable, indicates that
               real wkspace must be allocated in factor)

  UpperFill, LowerFill   class vars
*/

{
  int *lrow, *list;
  int i,   ii, iup, ilow;
  int first, iold, iend, num, next;
  int MAXINT;
  int nrow_tmp;
  int itemp, nrow;
  int *int_temp ;
 

  /* check validity of lorder, invord   */

  for (i=0; i<n;i++) {
    if (lorder[invord[i]] != i || porder[invpord[i]] != i) {
      throw General_Exception("invalid invord\n");
    }
  }

  /* allocate array of structures, each structue points to row */

  if (rowsp == NULL) {
    throw General_Exception("storage error in sfac2\n");
  }
  if (Dspace == allocated || Symfac == symfac_done || Numfac == numfac_done) {
    clean();
    init();
  }

  /* allocate temp wkspace */

  list = NULL;
  lrow = NULL;
  int_temp = NULL;

  try{// start try block
    list = new int [n];
    assert( list != NULL);
    lrow =  new int [n];
    assert( lrow != NULL);
    int_temp =  new int [n];
    assert( int_temp != NULL);
  }// end try block
  catch(std::bad_alloc){// catch
    delete [] list; list = NULL;
    delete [] lrow; lrow = NULL;
    delete [] int_temp; int_temp = NULL;
    throw;
  }// end catch

  try{ // start try block
    MAXINT = 2*n;
    *ier = 0;
    for (i=0; i < n; i++){
       list[i] = n+1;
       lrow[i] = MAXINT; 
    }
 
    itemp = 0;
    *nzero = 0;
    for (i = 0; i < n; i++) {
      iold = lorder[i];
 
      iend = 0;
      for (ii=ia[iold]; ii <= ia[iold+1]-1; ii++) {
        int_temp[iend] =  invpord[ ja[ii] ] ;
        iend++;
      }

      /* sort entries */
      num = ia[iold+1] - ia[iold];
      assert(num == iend);
      if (num > 1) {
        shell( int_temp , num);
      }

      first = int_temp[ 0 ];
      for (ii= 1; ii < num; ii++) {
        list[ int_temp[ii-1] ] = int_temp[ii];
      }
      list[int_temp[num-1]] = n+1;

      /* load level ptrs into lrow() array */
 
      for (ii=ia[iold]; ii < ia[iold+1]; ii++) {
        lrow[invpord[ja[ii]]] = 0;
      }
 
      merge2(i, list, first, lrow, level);
 
      /* count up nonzeros in this row  */
      nrow=0;
      next = first;
      while (next != n+1) {
        nrow++;
        next = list[ next ];
      }
      rowsp[i].nz =  nrow;
      *nzero = *nzero + nrow;

      /* allocate space for this row */
      if (rowsp[i].jaf != NULL) {
        throw General_Exception("error: jaf not deallocated in syfac\n");
      }
      rowsp[i].jaf = new int [nrow];
      assert( rowsp[i].jaf != NULL );
      if (rowsp[i].lev != NULL) {
        throw General_Exception("error: lev not deallocated in syfac\n");
      }
      rowsp[i].lev =  new int [nrow];
      assert(rowsp[i].lev != NULL);
 
      next = first;
      nrow_tmp=0;
      while (next != n+1) {
        rowsp[i].jaf[nrow_tmp] = next;
        rowsp[i].lev[nrow_tmp] = lrow[next];
        lrow[next] = MAXINT;
        if (next == i) {
          rowsp[i].diag = nrow_tmp;
        }
        next = list[next];
        nrow_tmp++;
      }
      assert(rowsp[i].jaf[rowsp[i].diag] == i);
      assert(nrow_tmp == nrow);
    }

    /* free level ptr space */
    for (i=0; i< n; i++) {
      delete [] rowsp[i].lev ;
      rowsp[i].lev = NULL;
    }

  }// end try block
  catch(...){ // start catch
    delete [] int_temp; int_temp = NULL;
    delete []  lrow;    lrow = NULL;
    delete []  list ;   list = NULL;
    throw;
  }// end catch
  
  /* free integer temp */
  delete [] int_temp;
  delete []  lrow;
  delete []  list ;

  //
  // indicate that symfac done
  Symfac = symfac_done;

  // must redo numfac now
  Numfac = numfac_not_done;

  // must redo wkspace allocation
  Dspace = unallocated;

  iup = 0;
  ilow = 0;
  for (i=0; i< n; i++) {
    ilow += rowsp[i].diag;
    iup += ( rowsp[i].nz - 1 - rowsp[i].diag );
  }
  UpperFill = iup;
  LowerFill = ilow;

  return;
}


/* factor */

void scaler_ILU::factor( const int ia[], const int ja[], const double a[]) 
 
/* numeric factor assuming sym fac done */

/*
   input

    ia, ja      stucture of a
    lorder[]    lorder[new_order] = old_order (class variable)
    invord      invord[old_order] = new_order (class variable)
                     row ordering
    porder[]
    invpord[]   porder[new_order] = old_order
                invpord[old_order] = new_order
                         col ordering
    n           number of unknowns  (class variable)
    a[]         real values of a
    *rowsp      contains info about symbolic factors
                of the ILU (class variable)
    Dspace      = allocated -> do not allocate real space
                for ILU, assume already done (i.e. prevous call
                to factor with same ILU structure, or prev
                call to facdrp
    Symfac      class variable

   output

    *rowsp      now contains real ILU as well as symbolic
                (ptr to global structure) (class variable)
                space is allocated in Dspace = unallocated
                on entry, otherwise, space overwritten
     Dspace     flagged allocated if done here
     Numfac     (class variable)
*/
{
  int i, iold, id, ii, idd, iii;
  int  nrow, *list, *ptr_jaf;
  double mult, *row; 
  SizeScalerFactor *ptr_af;


  /* allocate temp wkspace */

  if (Symfac != symfac_done) {
    throw General_Exception("Error: symbolic ILU not done before call to numeric ILU\n");
  }

  list = NULL;
  row = NULL;
     
  try{ // try block
    list = new int [n];
    assert( list != NULL);
    row = new double [n];
    assert( row != NULL);
  }// end try block
  catch(std::bad_alloc){// catch block
    delete [] list; list = NULL;
    delete [] row; row = NULL;
  }// end catch block

  for (i=0; i<n; i++) {
    row[i] = 0.0;
    list[i] = n+1;
  }

  try{ // start try block
 
    /* loop over rows */
    for (i=0; i<n; i++) {// start loop over rows

      iold = lorder[i];

      /* load orig matrix elements */
      for (ii=ia[iold]; ii<ia[iold+1]; ii++) {
        row[invpord[ja[ii]]] = a[ii]; // new col ordering
      }

      /* load markers */
      for (ii=0; ii<rowsp[i].nz; ii++) {
        list[ rowsp[i].jaf[ii]  ] = i;       
      }

      /* now eliminate */
      for (ii=0; ii<rowsp[i].diag; ii++) {

        id = rowsp[i].jaf[ii];
 
        /* get multiplier */
        mult = row[id] / rowsp[ id ].af[ rowsp[id].diag ];
        row[id] = mult ;
 
        ptr_jaf = rowsp[id].jaf;
        ptr_af = rowsp[id].af;
        for (iii=rowsp[id].diag+1; iii<rowsp[id].nz ; iii++) {
          idd = ptr_jaf[iii];
          if (list[idd ] == i) {
            row[idd ] -=  mult * ptr_af[ iii  ];
          }
        }
      }
 
      /* end of elimination for row */

      /* gather and reset row[] array */

      nrow = rowsp[i].nz;
      if (Dspace == unallocated) {
        if (rowsp[i].af != NULL) {
          throw General_Exception("error in scaler factor: af not deallocated\n");
        }
        rowsp[i].af = new SizeScalerFactor [nrow];
        assert(rowsp[i].af != NULL);
      }

      for (ii=0; ii<rowsp[i].nz ; ii++) {
        rowsp[i].af[ii] = row[rowsp[i].jaf[ii]];
        row[rowsp[i].jaf[ii]] = 0.0;
      }

      /* end of loop over rows */
 
    }// end loop over rows

  }// end try block
  catch(...){ // start catch block
    delete [] row; row = NULL;
    delete [] list; list = NULL;
    throw;
  }// end catch block

  /* deallocate temp wkspace */
  delete []  row ;
  delete []  list ;

  // reset Dspace to indicate that space has been
  // allocated
  if (Dspace == unallocated) {
    Dspace = allocated;
  }
 
  Numfac = numfac_done;
 
  return;
}


void scaler_ILU::solve(double x[], const double b[]) 

/*     forward and back solve  
         solve LU x = b

input:

   b[]       rhs
   n         number of unknowns (class variable)
   lorder[]  lorder[new_order] = old_order (class variable)
             row orderinge
   *rowsp    ptr to structure containing ILU
             (class variable)

   porder[new_order] = old_order   col ordering

output:

    x[]
*/
{
  int i, ii, *ptr_jaf;
  double *temp; 
  SizeScalerFactor *ptr_af;


  if (Numfac == numfac_not_done) {
    throw General_Exception("error: solve called with no factor\n");
  }
 
  /* forward solve:  Lz = b
     (L has unit diagonal)      */
 

  /* allocate temp */
  temp = new double [n];
  assert( temp != NULL);

  for (i=0; i<n; i++) {
    x[i] = b[ lorder[i] ]; // row ordering

    ptr_jaf = rowsp[i].jaf;
    ptr_af = rowsp[i].af;
    int upper = rowsp[i].diag;
    for (ii=0; ii< upper ; ii++) {
      x[i] -= ptr_af[ii] * x[ ptr_jaf[ii] ];
    }

   }
 
   /*  back solve Ux = z   */
   /*  (U does not have unit diag)    */
   for (i=n-1; i >= 0; i--) {

      ptr_jaf = rowsp[i].jaf;
      ptr_af = rowsp[i].af;

      int lower = rowsp[i].diag+1;
      int upper = rowsp[i].nz ;

      for (ii = lower; ii < upper ; ii++) {
        x[i] -=  ptr_af[ii] * x[ ptr_jaf[ii] ];
      }
      x[i] = x[i] /rowsp[i].af[ rowsp[i].diag ];
            
   }
 
   /* reorder */
   for (i=0; i< n; i++) {
     temp[ porder[i] ] = x[i]; // use col ordering
   }
 
   for (i=0; i<n; i++) {
     x[i] = temp[i];
   }
 
   delete [] temp;
;

   return;
}


/* free space for ilu in after finished  */

void scaler_ILU::clean(void)

/*
   input: (class variables) 
     n    number of unknowns
     *rowsp  structure containing ILU
*/
{
  delete [] rowsp;
  rowsp = NULL;
       
  return;
}

 
void scaler_ILU::merged(int i, int list[], int first, double arow[], 
double dropt, double odiag[])
/*
  input
   
     i     current row being factored (Dolittle form)
     list[]  implied linked list of nonzeros in current
             row.  On entry, only original nonzero columsn
             in this row.  1st nonzero column = fisrst
             , list[fisrt], list[ list[first]], ...., until
             list[  [....] ] = n+1 -> signals end of list
    n       number of nonzeros
             (class variable)
    first   (see list[]) fist nonzero in this row
    arow[]  original nonzeros (real numbers) for this
            row.  This is a full storage vector.  Nonzeros
            in original positions as given by list[].
            All other values = 0.
    dropt   drop tol

    odiag   sqrt of orig diag (DropType = 0)
            max orig row      (DropType = 1)
            orig diag         (DropType = 2)

    *rowsp  ptr to global structure containing info about
            previous rows, (real vals and nonzero structure)


    DropType   class var =0  sqrt( row * col )
                            = 1 max orig row
                            =2 orig diag
                             3 current diag
 
  output
 
    list[]   now contains all nonzeros in this factored row i
             (implied linked list)
    arow[]    numerical values of nonzero factors in this row
              arow = 0 at postions not given by list[]
*/
{
   int row, next, oldlst,  ii;
   int nxtlst;
   double mult, temp, test, drop_temp;


// drop test temps
                if( DropType == 0){
                       drop_temp = dropt *
                           odiag[i]  ;
                }
                else if( DropType == 1){
                       test  = dropt * odiag[i];
                }
                else if (DropType == 2){
                       test  = dropt * odiag[i];
                }


      next = first;

   while (  next < i) {

        oldlst = next;
        nxtlst = list[next];
        row = next;
 
        mult = arow[ oldlst] 
             / rowsp[row].af[ rowsp[row].diag ];
 
        arow[ oldlst ] = mult;

//      temp for drp test

                if (DropType == 3){
                     test  = dropt * fabs( arow[i]);
                }

/*          scan row "row" of U   */
 
    for( ii = rowsp[row].diag+1 ; ii <= rowsp[row].nz - 1 ; ii++) {
 
/*     skip thru linked list until we find possible insertion pt */

           while( rowsp[row].jaf[ii] > nxtlst) {
                 oldlst = nxtlst;
                 nxtlst = list[oldlst]; 
             }
                
           if( rowsp[row].jaf[ii] < nxtlst ){

                temp = mult  * rowsp[row].af[ii];
                if( DropType == 0){
                       test =  drop_temp *
                            odiag[ rowsp[row].jaf[ii] ] ;
                }
 
                if( fabs( temp) > test){
                  list[ oldlst] = rowsp[row].jaf[ii];
                  list[ rowsp[row].jaf[ii]  ] = nxtlst;
                  oldlst = rowsp[row].jaf[ii];
                  arow[ oldlst] -= temp;
                }
             }
           else {
              oldlst = nxtlst;

              arow[ oldlst ] -= mult*rowsp[row].af[ii];
              nxtlst = list[oldlst ];
            }
        }
           next = list[ next ];

  }
      return  ;
}


void scaler_ILU::facdrp(const int ia[], const int ja[], double a[],  
double dropt, int *nzero, int *ier)
/*
   input
   
      ia, ja, a    struct of a, and elements
      n            number of unkowns (class var)
      lorder       lorder(new_norder) = old_order
                      (class var)
      invord      invord[ lorder[i] ] = i  (class var)
                    (invord[ old_order] = new_order
      dropt       drop tol
      
      Drop[Type   class var =0  sqrt( row * col )
                            = 1 max orig row
                            =2 max orig col
                             3 current diag
  output

      *nzero      number of nonzeros in factors
      *ier        !=0 error in storage allocation
      *rowsp      class structure containg factors
      Dspace      class variable, which indicates
                  if real ILU space has been allocated
      Numfac, Symfac   flags to indicate the symbolic and
                       numeric factorization done
                   (ILU_class variables)
      UpperFill, LowerFill (class vars)
*/
{
      int  *list;
      int i,   ii;
     int  first, iold, iend, num,  next;
      int   nrow_tmp, ichek;
      int itemp, nrow;
      int *int_temp,  ilow, iup;
      double  *odiag, *arow, dmax;

      list = NULL; int_temp = NULL;
      odiag = NULL;  arow = NULL;
 

  try{ // start try block

/*     check validity of lorder, invord   */

      for(i=0; i<=n-1;i++){
         if( lorder[ invord[i] ] != i){
           throw General_Exception("invalid invord\n");
         }
      }

 
   if( rowsp == NULL){
        throw General_Exception("storage error in facdrp\n");
   }
    if( Dspace == allocated || Symfac == symfac_done ||
         Numfac == numfac_done){
       clean();
       init();
    }



/*     allocate temp wkspace      */

      list = new int [n];
         assert( list != NULL);
      arow = new double [n];
         assert( arow != NULL);
      odiag = new double [n];
         assert( odiag != NULL);
      int_temp = new int [n];
         assert( int_temp != NULL);



 
 
       *ier = 0;
 
 
       for (i=0; i <= n-1; i++){
          list[i] = n+1;
          arow[i] =  0.0;
       }
 
/*      find diag entries for odiag array   */


       for (i=0; i <= n-1; i++){
         ichek = 0;
         iold = lorder[i];
         dmax = 0.0;
         for(ii = ia[iold]; ii<= ia[iold+1]-1; ii++){
           dmax = max( dmax, fabs( a[ii] ) );
                if( ja[ii] == iold){
                   ichek = 1;
                   odiag[ i ] = fabs ( a[ii] );
                  }
         }
         assert( ichek == 1);
         if( DropType == 0){
            odiag[i] = sqrt(odiag[i]);
          }
         else if( DropType == 1){
            odiag[i] = dmax;
           }
       }
 
 
       itemp = 0;
       *nzero = 0;
     for( i = 0; i <= n-1; i++ ) {
         iold = lorder[i];
 
         iend = 0;
         for( ii=ia[iold]; ii <= ia[iold+1]-1; ii++){
               int_temp[iend] =  invord[ ja[ii] ] ;
               iend++;
         }

/*           sort entries */

         num = ia[iold+1] - ia[iold];
        
         assert( num == iend);
         if( num > 1){
             shell( int_temp , num);
         }

          first = int_temp[ 0 ];
          for( ii= 1; ii <= num-1; ii++){
                  list[ int_temp[ii-1] ] = int_temp[ii];
          }
          list[ int_temp[num-1]  ] = n+1;


/*       load orig entries into row  array */
 
          for( ii=ia[iold]; ii <= ia[iold+1]-1; ii++){
                arow[ invord[ ja[ii] ] ] =  a[ ii ];
          }
 
            merged( i,  list, first, arow, dropt, odiag);
 

 
/*                count up nonzeros in this row  */

           nrow=0;
           next = first;
           while(next != n+1){
              nrow++;
              next = list[ next ];
           }
            rowsp[i].nz =  nrow;
            *nzero = *nzero + nrow;

/*      allocate space for this row     */

            if( rowsp[i].jaf != NULL){
            throw General_Exception("error: jaf not dealocated in facdrp\n");
            }
            rowsp[i].jaf = new int [nrow];
               assert( rowsp[i].jaf != NULL );
            if( rowsp[i].af != NULL){
             throw General_Exception("error: af not dealocated in facdrp\n");
            }
            rowsp[i].af = new SizeScalerFactor [nrow];
               assert( rowsp[i].af != NULL);
              
 
           next = first;
           nrow_tmp=0;
           while( next != n+1){
                 rowsp[i].jaf[nrow_tmp] = next;
                 rowsp[i].af[nrow_tmp] = arow[next];
               arow[next] = 0.0;
               if( next == i){
                    rowsp[i].diag = nrow_tmp;
               }
               next = list[next];
               nrow_tmp++;
            }
               assert( rowsp[i].jaf[ rowsp[i].diag ] == i);
               assert( nrow_tmp == nrow);


    }


  }// end try block
  catch(...){// catch
      delete [] int_temp; int_temp = NULL;
      delete [] list; list = NULL;
      delete [] odiag; odiag = NULL;
      delete [] arow; arow = NULL;
      throw;
   }// end catch

/*       free integer temp   */

         delete [] int_temp;
         delete []  list ;

         delete []  odiag;
         delete []  arow ;
     
/*      count up number of entries in lower and uppper
       factors
*/

        iup = 0;
        ilow = 0;
      for(i=0; i<=n-1; i++){
        ilow += rowsp[i].diag;
        iup += ( rowsp[i].nz - 1 - rowsp[i].diag );
      }

//      printf("lower fill %d upper fill %d\n",
 //              ilow, iup);

        LowerFill = ilow;
        UpperFill = iup;
       
//    signal that real ILU space has been allocated
//    can be used to call factor without calling
//    the symbolic ILU, i.e. using drop ILU structure

     Dspace = allocated;

     Symfac = symfac_done;
     Numfac = numfac_done;

      return;
}
      

void scal(const int n, // number of unknowns
          int* ia, int* ja, 
          double* a, // orig matrix
          double* b, // right hand side
          double* scal_fac, // scaling factors
          int scaltype // = 0 diag
                       // = 1 max row
         )

//    scale matrix and rhs by inverse diag
//    alter matrix and rhs
//    store scaling factor in scal_fac
//    can be used for multiple rhsides
{ 
  const double eps = 1.0e-300;
  int i;

  for (i=0; i<=n-1; i++) { // start n-loop
    int j;
    int found = -1;
    double dmax = 0.0;
    for (j=ia[i]; j<= ia[i+1]-1; j++) {// start j-loop

      if (scaltype == 0) {
        if (ja[j] == i) {
          found = j;
          break;
        }
      }
      else{
        // scal by max diag
        dmax = max( fabs(a[j]), dmax);
      }
    }// end j-loop

    if (scaltype == 0) { // no pivoting
      assert(found != -1);
      scal_fac[i] = 1.0/( a[ found] + eps);
    }
    else {
      scal_fac[i] = 1.0/(dmax + eps);
    }

    for (j=ia[i]; j<= ia[i+1]-1; j++) { 
      a[j] *= scal_fac[i];
    }

    b[i] *= scal_fac[i];
          
  }// end n-loop
}// end function scal


base_ILU::base_ILU(int n_in, int lord_in[], int* pord_in)
{
  lorder = NULL; invord = NULL;
  porder = NULL; invpord = NULL;

  int i;
  n = n_in;
  lorder = new int [n];
  assert( lorder != NULL);
  invord = new int [n];
  assert( invord != NULL);

  for (i=0;i<=n-1;i++) {
    lorder[i] = lord_in[i];
    assert( lorder[i] <n && lorder[i] >= 0);
  }
  for (i=0;i<=n-1;i++) {
    invord[ lorder[i] ] = i;
  }

  porder = new int [n];
  assert(porder != NULL);
  invpord = new int [n];
  assert(invpord != NULL);


  for (i=0; i<n; i++) {

    if (pord_in != NULL) {
      porder[i] = pord_in[i];
    }
    else {
      porder[i] = lorder[i];
    }
    invpord[ porder[i] ] = i;
  }

  Numfac = numfac_not_done;
  Symfac = symfac_not_done;
  Dspace = unallocated;
  UpperFill = 0;
  LowerFill = 0;
} // end base_ILU::base_ILU( int n_in, int lord_in[])


void scaler_ILU::init(void)
{
  rowsp = NULL;

  rowsp = new  scaler_sparse_row [ n ] ;
  assert( rowsp != NULL);
  Numfac = numfac_not_done;
  Symfac = symfac_not_done;
  Dspace = unallocated;
  DropType = 0;
  PivType = 0;

}// end scaler_ILU::init( void)


void scaler_ILU::facdrp2(const int ia[], const int ja[], double a[],  
double dropt, int *nzero, int *ier)
/*
   input
   
      ia, ja, a    struct of a, and elements
      n            number of unkowns (class var)
    row ordering
      lorder       lorder(new_norder) = old_order
                      (class var)
      invord      invord[ lorder[i] ] = i  (class var)
                    (invord[ old_order] = new_order

     col ordering
                    porder[new_order] = old_order
                    invpord[old_order] = new_order

      dropt       drop tol
      
      int Drop[Type   class var =0  sqrt( row * col )
                            = 1 max orig row
                            =2 max orig col
                             3 current diag
      int PivType  = 0  // no pivoting
                       // = 1 pivot
  output

      *nzero      number of nonzeros in factors
      *ier        !=0 error in storage allocation
      *rowsp      class structure containg factors
      Dspace      class variable, which indicates
                  if real ILU space has been allocated
      Numfac, Symfac   flags to indicate the symbolic and
                       numeric factorization done
                   (ILU_class variables)
      UpperFill, LowerFill (class vars)
*/
{
/*     check validity of lorder, invord   */

   { int i; // start block
      for(i=0; i<=n-1;i++){
         if( lorder[ invord[i] ] != i ||
             porder[ invpord[i] ] != i ){
             throw General_Exception("invalid invord\n");
         }
      }
   } // end block

 
{ //check block
   if( rowsp == NULL){
         throw General_Exception("storage error in facdrp\n");
   }
    if( Dspace == allocated || Symfac == symfac_done ||
         Numfac == numfac_done){
       clean();
       init();
    }

 }// end check block



/*     allocate temp wkspace      */

      double* arow = NULL;
      double* odiag = NULL;
      int* marker = NULL;
      int* lower = NULL;
      int* upper = NULL;

 try{ // start try block
      
       arow = new double [n];
         assert( arow != NULL);
       odiag = new double [n];
         assert( odiag != NULL);

       marker = new int [n];
          assert( marker != NULL);
       lower = new int [n];
          assert( lower != NULL);
       upper = new int [n];
          assert(upper != NULL);


// 

 
 
       *ier = 0;
 
 
// initialize full length arrays

    { int i;
       for (i=0; i <= n-1; i++){
          arow[i] =  0.0;
          marker[i] = -1;
          odiag[i] = 0.0;
     }
   }
 

// find data for use in drop test


    { // start block
         int i;
       for (i=0; i <= n-1; i++){// loop over i
         int ichek = 0;
         int iold = lorder[i]; // row pivot
         double dmax = 0.0;
         int ii;
         for(ii = ia[iold]; ii<= ia[iold+1]-1; ii++){ // loop over
                                              // neighbors
           dmax = max( dmax, fabs( a[ii] ) );
           if( DropType == 0){ // diagoanl
                if( ja[ii] == iold){
                   ichek = 1;
                   odiag[ i ] = fabs ( a[ii] );
                  }
           }// end diago

          if( DropType == 2){// max col
                 odiag[ ja[ii] ] = max( odiag[ ja[ii] ],
                                        fabs( a[ii] ));
            }// end max col

         }  // end loop over neighs

       
        if( PivType == 0){ 

         if( DropType == 0){
            assert( ichek == 1);
            odiag[i] = sqrt(odiag[i]);
          }
         else if( DropType == 1){
            odiag[i] = dmax;
           }

         }
         else{
            odiag[i] = dmax;
            DropType = 1;
          }

       }// end i loop

     } // end block
 
 
  {// main loop block
       *nzero = 0;
       int i; // row index, new ordering

     for( i = 0; i <= n-1; i++ ) { // main loop
         int iold = lorder[i];
 
      int count_up = -1;
      int count_low = -1;

      { int ii; // load orig entries into row

         for( ii=ia[iold]; ii <= ia[iold+1]-1; ii++){// loop over neighs

             int id = invpord[ ja[ii] ]; // new col ordering

              arow[ id ] =  a[ ii ];
             if( id < i ){
                count_low++;
                lower[count_low] = id;
             }
             else{
               count_up++;
               upper[count_up] = id;
             }
             marker[ id ] = i;
         }// end loop over neighs

       } // end load orig entries

// on entry. prevous rows are in the new row order
//  (rowsp[k], but cols are:
//        L (j < diag) rowsp[k].jaf[j] in New col order
//        U (j >= diag) rowsp[k].jaf[j] in OLD col order
//
//   arow input/output in NEW col order (but this is only
//     the current new col order, which may change, sow
//     we are going to store this in the old orig col order

//  lower[], upper[] hold col indicies for this row
//   New order
// on exit from elimrow, lower[] is sorted, upper[] is not sorted

     elimrow(i, odiag, arow, marker, count_up, count_low,
             lower, upper, dropt);
 
/*                count up nonzeros in this row  */

//     chek to make sure that there is a diag
//      entry

           assert( count_up >= 0);

           int nrow= (count_up +1) + (count_low+1);
            rowsp[i].nz =  nrow;
            rowsp[i].diag = count_low+1;

            *nzero = *nzero + nrow;

/*      allocate space for this row     */

            if( rowsp[i].jaf != NULL){
             throw General_Exception("error: jaf not dealocated in facdrp\n");
            }
            rowsp[i].jaf = new int [nrow];
               assert( rowsp[i].jaf != NULL );
            if( rowsp[i].af != NULL){
             throw General_Exception("error: af not dealocated in facdrp\n");
            }
            rowsp[i].af = new SizeScalerFactor [nrow];
               assert( rowsp[i].af != NULL);
              
       // load in new entries in NEW col order

            {// load block
               int k;
               int count_tot= -1;

               for(k=0; k<=count_low; k++){
                  count_tot++;
                  rowsp[i].af[count_tot] = arow[ lower[k] ];
                        arow[ lower[k] ] = 0.0;
                  rowsp[i].jaf[count_tot] = lower[k];
               }

               int isave = -1;


         // upper[] is in the new ordering, but is
         // unsorted
         // Note: ordering for lower[] cannot change
         // but ordering for upper[] may change due
         // to pivoting in future rows

               for(k=0; k<=count_up; k++){
                  count_tot++;
                  if( upper[k] == i){
                        isave = count_tot;
                   }
                  rowsp[i].af[count_tot] = arow[ upper[k] ];
                        arow[ upper[k] ] = 0.0;
                  rowsp[i].jaf[count_tot] = upper[k];

               }

                assert( isave != -1);

              // swap to diagonal
              // so that row[i].jaf[ diag ] = i
              // this col order cannot be changed now,
              // so this must be the diagonal in new ordering

              {   int diag = rowsp[i].diag;
                 int itemp = rowsp[i].jaf[diag];
                 rowsp[i].jaf[diag] = rowsp[i].jaf[isave];
                 rowsp[i].jaf[isave] = itemp;
              
                 double temp = rowsp[i].af[diag];
                   rowsp[i].af[diag] = rowsp[i].af[isave];
                   rowsp[i].af[isave] = temp;
              }

               assert( (count_tot+1) ==  nrow);

     // NOW, convert U back to original col ordering
     // this is because the ordering may change if future
     // rows due to pivoting, so we save back in orig ordering

              { // relabel upper part of rowsp.jaf[]
                int k;
                for(k=rowsp[i].diag; k<rowsp[i].nz; k++){
                   int old_order = porder[ rowsp[i].jaf[k] ];
                     rowsp[i].jaf[k] = old_order;
                 }
              }

                

            }// end load block


    } // end main loop
  }// end main loop block


   { // fix up rowsp[i].jaf so that is in new ordering
     // U only
         int i;
         for(i=0; i<n; i++){
             int k;
                for(k=rowsp[i].diag; k<rowsp[i].nz; k++){
                   int new_order = invpord[ rowsp[i].jaf[k] ];
                     rowsp[i].jaf[k] = new_order;
                 }
          }
   }// end fix


}// end try block
 catch(...){
       delete []  odiag; odiag = NULL;
       delete []  arow ; arow = NULL;
       delete []  marker ; marker = NULL;
       delete []  lower ;  lower = NULL;
       delete [] upper ;   upper = NULL;
       throw;
   }// end catch


/*       free integer temp   */


         delete []  odiag;
         delete []  arow ;

       delete []  marker ;
       delete []  lower ;
       delete [] upper ;

     
/*      count up number of entries in lower and uppper
       factors
*/

    { //stat block
        int iup = 0;
        int ilow = 0;
        int i;
      for(i=0; i<=n-1; i++){// i-loop
        ilow += rowsp[i].diag;
        iup += ( rowsp[i].nz - 1 - rowsp[i].diag );
      } // end i-loop

//      printf("lower fill %d upper fill %d\n",
 //              ilow, iup);

        LowerFill = ilow;
        UpperFill = iup;
       
//    signal that real ILU space has been allocated
//    can be used to call factor without calling
//    the symbolic ILU, i.e. using drop ILU structure

     Dspace = allocated;

     Symfac = symfac_done;
     Numfac = numfac_done;
    } //end stat block



      return;
}// end facdrp2
      

void scaler_ILU::elimrow( const int i , // row being elimiated
                          double* odiag, // dtop tol test data
                          double* arow, // full length arrey of
                                       // nonzeros for this row
                          int* marker, // nonzero in col_j ->
                                       // marker[j] = i
                          int& count_up,  
                          int& count_low, // col entries for this row
                          int* lower,     // for col < diag,
                          int* upper,      // lower[0,..,count_low]
                                          // for col >=diag
                                          // upper[0,..,count_up]
                                          // on entry, unordered
                                          // on exit, lower[] sorted
                                          // upper[] unsorted
                          double dropt // drop tol
                        )


//
//  int PivType  = 0  // no pivoting class var
//                       // = 1 pivot

 {
          double drop_temp;
          double test;

               if( DropType == 0){
                       drop_temp = dropt *
                           odiag[i]  ;
                }
                else if( DropType == 1){
                       test  = dropt * odiag[i];
                }
                else if (DropType == 2){
                       test  = dropt * odiag[i];
                }

//Note:   upper U for previous rows (rowsp[j].jaf[ ]
//       is in ORIGINAL col ordering
//       
//     however, arow[], upper[], lower[] are in the
//     NEW col ordering

   int k; 

   for(k=0; k<=count_low; k++){ // eliminate lower factor 
                                // entries
     




    { // start block
      int ipiv = n+1;
      int j; 
       int isave = -1;
        for(j=k; j<= count_low; j++){ // find smallest unelim col index
                                      // entry
          if( lower[j] < ipiv){
            ipiv = lower[j];
            isave = j;
           }
         } // end find smallest unelim col index

       assert( isave != -1);

      // swap entry
      { int itemp = lower[isave];
        lower[isave] = lower[k];
        lower[k] = itemp;
      }

    } // end block


      {// elimination block
           int row = lower[k];

        double mult = 
           arow[ row] /rowsp[ row ].af[ rowsp[row].diag];
           arow[ row ] = mult;


        if (DropType == 3){
                 test  = dropt * fabs( arow[i]);
             }// current diag test

         

         int ii;

         for(ii = rowsp[ row].diag+1; 
                  ii <rowsp[row].nz; ii++){ // row  loop

                 int ifill = invpord[ rowsp[row].jaf[ii]]; 
                               // convert this to new col order

                 if( marker[ifill] == i){// nonzero entry here

                     arow[ifill] -= mult*rowsp[row].af[ii];
                    
                }// end nonzero here already

                else{// new fill-in

                   double temp =  mult  * rowsp[row].af[ii];
                   if( DropType == 0){
                       test =  drop_temp *
                            odiag[ invord[ rowsp[row].jaf[ii] ]  ] ;
                                   // convert to new row order

                     }
                     else if( DropType == 2){
                           test = dropt*odiag[ rowsp[row].jaf[ii] ];
                         }
                    
                    if( fabs(temp) > test ){ // insert new entry
                        
                            arow[ ifill ] -= temp;
                             if( ifill < i ){
                                count_low++;
                                lower[count_low] = ifill;
                               }
                               else{
                                   count_up++;
                                   upper[count_up] = ifill;
                                 }
                                 marker[ ifill ] = i;
                   }// end insert new entry


                 }// end new-fill in

         } // end row loop

                
    } // end elmination block


  } // end eliminate lower factor entries


   { // pivot block
      if( PivType == 0){ 
          return;
      }

      int isave_diag = -1;
      int isave_max = -1;

      int k;
      double dmax = 0.;
      for(k=0; k<=count_up; k++){ // loop over U

         if( fabs( arow[upper [k] ] ) > dmax){
              isave_max = k;
              dmax = fabs( arow[upper [k] ] );
          }
          if( upper[k] == i){
              isave_diag = k;
           }
       } // end loop over U

           // lower[], upper[] in new col order
           // upper[] is unsorted
  
             if( isave_max == -1){
                throw General_Exception("ILU is singular\n");
              }


       int iswap = upper[isave_max];

       {// swap row entry
          double temp = arow[i];
           arow[i] = arow[iswap];
           arow[iswap] = temp;
       }

       { // swap perm vectors
           int itemp = porder[i];
              porder[i] = porder[iswap];
              porder[iswap] = itemp;
 
              invpord[ porder[i] ] = i;
              invpord[ porder[iswap] ] = iswap;
       }

       

       if( isave_diag == -1){
         // no col entry on diagaonal
         //
         // simply relabel
            upper[isave_max] = i;
         }

      // if isave_diag != -1, then there already
      // is a diag and iswap entry, simply do nothing
      // to upper[]

     
   }// end pivot block


 
}// scaler_ILU:: end elimrow

}// end namespace SparseItObj{

