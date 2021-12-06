
#ifndef DEF_INC
#define DEF_INC


//  remove these for other compilers
#undef UNIX_BUILD_
#define UNIX_BUILD_
#undef OLD_IOSTREAM

// Set a flag to indicate whether we are using the Microsoft
// Visual C++ compiler.

//#define OLD_IOSTREAM
//#undef OLD_IOSTREAM

#ifndef UNIX_BUILD_
    #define VISUAL
    #define OLD_IOSTREAM
#endif

#ifdef UNIX_BUILD_
    #undef VISUAL
    #ifndef OLD_IOSTREAM
       #define NEW_HEADERS // new headers, unix build only
    #endif
#endif

// new/old iostream library
//



//   remove comments to kill assert cheks
//#define NDEBUG

#endif

