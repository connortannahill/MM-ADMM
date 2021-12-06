#include "def_compiler.h"
#include "Standard.h"
#include "SparseItUtil.h"


using namespace SparseItObj;
using namespace std;

/* prototypes   */

namespace SparseItObj{ 

    void rcm( const int ia[], const int ja[], int n, int lorder []);

    void levels( int ia[], int ja[], int n, int iroot, 
       int *nlvl, int ial[], int jal[] );

    void rcm1( int ia[], int ja [], int n, int iroot ,
        int nlvl, int ial[], int jal[], int lorder[]);
 
    void pseudo( int ia[], int ja[], int n, int ial[],
        int jal[], int* iroot, int* nlvl);

   void rcm1( int ia[], int ja[], int n, int iroot,
     int nlvl, int ial[], int jal[], int lorder[]);

  void symmetrize( const int*  ia, // input 
                   const int* ja , // input
                   const int n, 
                   int* & ia_new, 
                   int* & ja_new ); // output symmetrized ia-ja
/* function defs  */

void rcm( const int ia[], const int ja[], int n, int lorder[])
   {

/*    determine level structure of ia-ja graph

  input:
       ia, ja  structure of graph
        n       number of nodes
        nja     size of ja() array

NOTE:   stucture of matrix is symmetrized before
        ordering, does not alter original ia-ja


   output:

     lorder( new_order ) = old_order


*/

     int  i, id,  ii,  ichek,
          jj,  nlvl,  iroot;

      
/*     check that incidence matrix is symmetric */


      for(i = 0; i<= n-1; i++){
         for(ii = ia[i]; ii <= ia[ i+1]-1; ii++){
            id = ja[ii];
            ichek = 0;
            for( jj = ia[id]; jj <= ia[id+1]-1; jj++){ //jj-loop
               if( ja[jj] == i){
                  if( ichek != 0){
        throw General_Exception("error rcm: same index appears twice\n");
                   }
                   else{
                    ichek = 1;
                   }
               }
            } // end jj-loop
      //      if( ichek == 0){
      //        cout << "nonsymmetric structure\n";
      //        status = 1;
      //       exit( status ); 
      //    }

         } // end ii loop 
      } // end i-loopo

 { 
   int* ia_new = NULL; ; // temp symmetrized ia
   int* ja_new  = NULL;// temp symmetrized ja

    int* ial = NULL;
    int* jal = NULL;

  try{ // start try block

    symmetrize( ia, ja, n, ia_new, ja_new);

    int nja; // more temps 

      ial = new int [n+1];
        assert ( ial != NULL);
      nja =  ia_new[n] ;

      jal = new int [nja];
        assert( jal != NULL);


   iroot = 0;
   pseudo (ia_new , ja_new, n , ial, jal, 
           &iroot, &nlvl );
   rcm1( ia_new, ja_new ,  n,  iroot ,
         nlvl,  ial,  jal,  lorder);

   }// end try block
   catch(...){
       delete [] ia_new; ia_new = NULL;
       delete [] ja_new; ja_new = NULL;
       delete [] ial; ial = NULL;
       delete [] jal; jal = NULL;
       throw;
     }
   
    delete [] ia_new;
    delete [] ja_new; // delete temporaries

      delete[]  ial ;
      delete []  jal ;


 }// end block

      
  // chek valid ordering
  { int* ichek;
     ichek = new int [n];
       assert( ichek != NULL);
      int ii;
      for(ii=0; ii<n; ii++){ // start loop over nodes
         ichek[ii] = -1;
       }// end loop over nodes
      for(ii=0; ii<n; ii++){ // start loop over nodes
         int jj = lorder[ii];

//  debug***********
//       printf( " ii lorder[ii]: % d % d \n", ii, lorder[ii]);

         if( jj < 0 || jj >= n){
           throw General_Exception("error_1 in rcm ordering\n");
          } 

         if( ichek[jj] != -1){
            throw General_Exception("error_2 in rcm ordering\n"); 
          }
           ichek[jj] = +1;

       }// end loop over nodes

    delete [] ichek;
  }

   return;

  }

    void levels( int ia[], int ja[], int n, int iroot,
       int *nlvl, int ial[], int jal[] )
  {

/*    get level structure starting at node "iroot" 

   input:

        ia, ja  structure of graph
        n       number of nodes
        nja     size of ja() array
        iroot   starting node

     output:


         nlvl  number of levels in level struc
               starting at iroot
         ial(), jal()   ptrs which describe nodes
                in level structure
     

*/

    
    int *proces ,i, id, idd,  icount, nlold, jj;

    proces = new int [n];
      assert( proces != NULL);

/*   mark all nodes as unprocessed   */

    for( i=0; i <= n-1 ; i++){
       proces[i] = 1;
    }

    ial[0] = 0;
    icount = 0;
    jal[ icount ] = iroot;
    proces[ iroot ] = 0;
    nlold = 0 ;
    ial[ nlold + 1] = icount + 1;
    *nlvl = nlold;

    while( icount < n-1 ){
      *nlvl = *nlvl + 1;

       if( *nlvl >= n){
          throw General_Exception(" rcm ordering will not work: matrix reducible\n");
        }

      for( i = ial[nlold]; i <= ial[nlold+1]-1; i++){
         id = jal[i];
         for(jj = ia[id]; jj <= ia[id+1]-1; jj++){
            idd = ja[jj];
            if( proces[idd] != 0){
               icount++;
               jal[icount] = idd;
               proces[ idd ] = 0;
            }
         }
      }
      ial[ *nlvl +1] = icount + 1;
      nlold = *nlvl;
    }

    assert( 
        (ial[ *nlvl + 1] - 1) == (n-1) 
          );
 
    delete []  proces ;

    return;
  }

 void pseudo( int ia[], int ja[], int n, int ial[],
        int jal[], int* iroot, int* nlvl)

 {

/*     get pseudo periperhal node   

   input:

        ia, ja  structure of graph
        n       number of nodes
        nja     size of ja() array
        iroot   starting node

   output:

      iroot   periperal node found
      nlvl    number of levels in level struc
      ial, jal  nodes in level structure
*/

   int nlold, icount, id, mindeg, i;

   levels(  ia,  ja,  n, *iroot, 
       nlvl,  ial,  jal );

   if( *nlvl == 1 || *nlvl == n){
      return;
   }

   nlold = -1;
   while ( *nlvl > nlold && *nlvl < n-1){
     nlold = *nlvl;
     mindeg = n+1;
     *iroot = jal[ ial[ *nlvl] ];
     for( i = ial[ *nlvl ]; i <= ial[ *nlvl + 1]-1; i++){
      id = jal[i];
      icount = ia[id+1] - ia[id];
      if( icount < mindeg) {
        mindeg = icount;
        *iroot = id;
      }
     }
       assert( mindeg != n+1);
 

     levels(  ia,  ja,  n, *iroot, 
       nlvl,  ial,  jal );


   }


   return;
 }


 
   void rcm1( int ia[], int ja[], int n, int iroot,
     int nlvl, int ial[], int jal[], int lorder[])
{ 
/*
   input:

        ia, ja  structure of graph
        ial, jal   level structure
                   starting at iroot
        n       number of nodes
        nja     size of ja() array
        iroot   starting node

     output:

      lorder(new_order) = old_order

*/

    int *proces, deg, i, ii, jj, id, ilvl, icold;
    int  icount, idd, itemp, iii;

/*    temp vector   */

    proces = new int [n];
     assert( proces != NULL);

    for (i=0 ; i<= n-1; i++){
      proces[i] = 1;
      lorder[i] = -1;
    }

    proces[ iroot ] = 0;
    icount = 0;
    lorder[ icount ] = iroot;
    
    for( ilvl=0 ; ilvl <= nlvl; ilvl++){


       for( iii = ial[ilvl] ; iii <= ial[ilvl+1]-1; iii++){
          icold = icount;
          id = lorder[ iii ];
          assert( id != -1);
          for( jj=ia[id]; jj<= ia[id+1]-1; jj++){
             idd = ja[jj];
             if( proces[idd] != 0){
                icount++;
                assert( lorder[icount] == -1);
                lorder[icount] = idd;
                proces[ idd ] = 0;
             }
          }
/*   sort nodes in order of incresing degree */

          for( ii=icold+2; ii <= icount; ii++){
             deg = ia[ lorder[ii] + 1 ] - ia[ lorder[ii] ];
             for ( jj= ii-1; jj >= icold+1; jj--){
                 if( deg >= ( ia[ lorder[jj]+1 ] - ia[ lorder[jj] ] ) ){
                      break;
                 }
                 itemp = lorder[jj];
                 lorder[jj] = lorder[jj+1];
                 lorder[jj+1] = itemp;
             }
          }

/*          end of sort   */

/*           check ordering    */

           for(ii=icold+2; ii <= icount; ii++){
             assert( 
                  ( ia[ lorder[ii] + 1 ] - ia[ lorder[ii] ])
                     >= (ia[ lorder[ii-1] + 1 ] - ia[ lorder[ii-1] ])
                    );
           }
        }
    }

/*   chek total  */

    assert( icount == n-1);
 
/*     reverse order, use ial as temp   */

    for(i=0, itemp=n-1; i<=n-1; i++, itemp--){
      ial[ itemp ] = lorder[ i ];
    }

    for(i=0; i<=n-1; i++){
         lorder[ i ] = ial[i];
    }
 
/*    check ordering     */
 
    for(i=0; i<=n-1; i++){
      ial[i] = -1;
    }
    for(i=0; i<=n-1; i++){
       assert( lorder[i] <= n-1 && lorder[i] >= 0);
       assert( ial[ lorder[i] ] == -1);
       ial[ lorder[i] ] = +1;
    }

    delete [] proces;

 return;
}

void symmetrize( const int* ia, // input 
                   const int* ja , // input
                   const int n, 
                   int* & ia_new, 
                   int* & ja_new )// output symmetrized ia-ja
 {

  // if necessary, ensure that ia-ja new rep a symmetric
  // incidence matrix

      int* count = new int [n];
        assert( count != NULL);
      
    { int i;
       for(i=0; i<n; i++){ count[i] = 0; }
    }
 

      {int i,j; // start block
        for(i=0; i<n; i++){// i loop
          for(j=ia[i]; j<ia[i+1]; j++){ // j-loop
              count[i]++;

              int id = ja[j];
                {// start block
                  int ichek = -1;
                  int k;
                   for(k=ia[id]; k< ia[id+1]; k++){ // k-loop
                      if( ja[k] == i){
                            ichek = 1; // found neighbour
                            break;
                       }
                    }// end k-loop

                   if( ichek == -1 ){ // no neighbour, leave space
                      count[id]++;
                   }
                }// end block

            } // end j loop
        }// end i loop
      }// end  block

      ia_new = new int[n+1];
        assert( ia_new != NULL);
 
      {int i; // start block

         ia_new[0] = 0;
         for(i=0; i<n; i++){
            ia_new[i+1] = ia_new[i] + count[i];
          }// end loop

      }// end block


     ja_new = new int [ ia_new[n] ];
         assert( ja_new != NULL);

     {int i;
       for(i=0; i<n; i++){
           count[i] = 0;
        }
     }
     

      {int i,j; // start block

        for(i=0; i<n; i++){// i loop
          for(j=ia[i]; j<ia[i+1]; j++){ // j-loop
              if( ia_new[i]+ count[i] > ia_new[i+1]-1 ){
                     throw General_Exception(" error_3 in rcm\n");
                }

              ja_new[ia_new[i] + count[i] ] = ja[j];
                   count[i]++;

              int id = ja[j];
                {// start block
                  int ichek = -1;
                  int k;
                   for(k=ia[id]; k< ia[id+1]; k++){ // k-loop
                      if( ja[k] == i){
                            ichek = 1; // found neighbour
                            break;
                       }
                    }// end k-loop

                   if( ichek == -1 ){ // no neighbour, add it

                     if( ia_new[id]+ count[id] > ia_new[id+1]-1 ){
                       throw General_Exception(" error_4 in rcm\n");
                      }

                      ja_new[ ia_new[id] + count[id] ] = i;
                      count[id]++;
                   } // end add neighbour

                }// end block

            } // end j loop
        }// end i loop
      }// end  block

      
      delete [] count; // clean store
 }// end function symmetrize

 } // end namespace SparseItObj{ 
