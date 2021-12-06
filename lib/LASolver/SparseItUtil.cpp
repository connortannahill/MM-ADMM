
#include "def_compiler.h"



#include "Standard.h"


/*    prototypes  */

#include "SparseItUtil.h"

using namespace SparseItObj;
using namespace std;

 


 namespace SparseItObj{ 

void shell( int v[], int n)

/*    shellsort:  sort v[0],... v[*n_ptr-1]
      into increasing order */

{  int gap, i, j, temp;

    for (gap = n/2; gap > 0; gap = gap/2){
       for (i = gap; i < n; i++){
          for( j = i - gap; j>=0 && v[j]>v[j+gap]; j=j-gap){
              temp = v[j];
              v[j] = v[j+gap];
              v[j+gap] = temp;
          }
       }
    }
}


 int bsearchIter( const int* array, // array to be searched
                         const int array_size, // size of array
                         const int key // key to be searched for
                        )
// assumes array ordered in increasing size
 {
     int lo = 0;
     int mid;
     int hi = array_size -1;

     while( lo <= hi){
       mid = (lo + hi)/2 ;

       if( key < array[ mid ]){ // in lower half
         hi = mid -1;
        } // end lower half

        else { // not lower half

             if( array[mid] < key){ // upper half
                 lo = mid + 1;
               }// end upper half
               else{ // found
                 return mid;
                }// end found

         }// not lower half

     }// endwhile

      return -1 ; // not in list

 }// end bsearchIter





#ifdef VISUAL
  int SparseItObj::freeStoreExhaust(size_t size ){
    throw std::bad_alloc();
 }

#else

void freeStoreExhaust(){
  cerr << " **************** NO MORE STORAGE *********" << endl;
  cerr <<" ***************** PROGRAM EXITING *********" << endl;
  cerr << "*******************************************" << endl;
  throw std::bad_alloc();
 }

#endif

 } // end namespace SparseItObj{
