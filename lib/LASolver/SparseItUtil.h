
#ifndef UTIL_INC
#define UTIL_INC

namespace SparseItObj {// start  namespace SparseItObj

class General_Exception {
public:
   const char* p;
   General_Exception(const char* q) {p=q;}
};


// utility functions, ie, max. min sort, etc

inline int min(const int i, const int j)
{
  if (i < j) {
    return i;
  }
  else {
    return j;
  }
} // end min()

inline int max(const int i , const int j)
{
  if (i > j) {
    return i;
  }
  else {
    return j;
  }
} // end max

void shell(int v[], int n); // shell sort in increasing order

inline double max(const double x, const double y)
{
  if (x > y)
    {return x;}
  else
    {return y;}
} // end max

inline double min(const double x, const double y)
{
  if (x < y)
    {return x;}
  else
    {return y;}
} // end min

int bsearchIter(const int* array,     // array to be searched
                const int array_size, // size of array
                const int key         // key to be searched for
               );
// assumes array ordered in increasing size

#ifdef VISUAL
  int freeStoreExhaust(size_t size);
#else
  void freeStoreExhaust();
#endif

}// end  namespace SparseItObj

#endif
