/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef __libmatlabTest_h
#define __libmatlabTest_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_libmatlabTest
#define PUBLIC_libmatlabTest_C_API __global
#else
#define PUBLIC_libmatlabTest_C_API /* No import statement needed. */
#endif

#define LIB_libmatlabTest_C_API PUBLIC_libmatlabTest_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_libmatlabTest
#define PUBLIC_libmatlabTest_C_API __declspec(dllexport)
#else
#define PUBLIC_libmatlabTest_C_API __declspec(dllimport)
#endif

#define LIB_libmatlabTest_C_API PUBLIC_libmatlabTest_C_API


#else

#define LIB_libmatlabTest_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libmatlabTest_C_API 
#define LIB_libmatlabTest_C_API /* No special import/export declaration */
#endif

extern LIB_libmatlabTest_C_API 
bool MW_CALL_CONV libmatlabTestInitializeWithHandlers(mclOutputHandlerFcn error_handler,
                                                      mclOutputHandlerFcn print_handler);

extern LIB_libmatlabTest_C_API 
bool MW_CALL_CONV libmatlabTestInitialize(void);

extern LIB_libmatlabTest_C_API 
void MW_CALL_CONV libmatlabTestTerminate(void);



extern LIB_libmatlabTest_C_API 
void MW_CALL_CONV libmatlabTestPrintStackTrace(void);


extern LIB_libmatlabTest_C_API 
bool MW_CALL_CONV mlxMatlabTest(int nlhs, mxArray *plhs[],
                                int nrhs, mxArray *prhs[]);


extern LIB_libmatlabTest_C_API bool MW_CALL_CONV mlfMatlabTest(int nargout
                                                               , mxArray** vec_out
                                                               , mxArray* vec_in);

#ifdef __cplusplus
}
#endif

#endif

