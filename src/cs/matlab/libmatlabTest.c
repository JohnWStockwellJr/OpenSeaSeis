/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <stdio.h>
#define EXPORTING_libmatlabTest 1
#include "libmatlabTest.h"
#ifdef __cplusplus
extern "C" {
#endif

extern mclComponentData __MCC_libmatlabTest_component_data;

#ifdef __cplusplus
}
#endif


static HMCRINSTANCE _mcr_inst = NULL;


#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultPrintHandler(const char *s)
{
  return mclWrite(1 /* stdout */, s, sizeof(char)*strlen(s));
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultErrorHandler(const char *s)
{
  int written = 0;
  size_t len = 0;
  len = strlen(s);
  written = mclWrite(2 /* stderr */, s, sizeof(char)*len);
  if (len > 0 && s[ len-1 ] != '\n')
    written += mclWrite(2 /* stderr */, "\n", sizeof(char));
  return written;
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libmatlabTest_C_API 
#define LIB_libmatlabTest_C_API /* No special import/export declaration */
#endif

LIB_libmatlabTest_C_API 
bool MW_CALL_CONV libmatlabTestInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler
)
{
  if (_mcr_inst != NULL)
    return true;
  if (!mclmcrInitialize())
    return false;
  if (!mclInitializeComponentInstanceWithEmbeddedCTF(&_mcr_inst,
                                                     &__MCC_libmatlabTest_component_data,
                                                     true, NoObjectType,
                                                     LibTarget, error_handler,
                                                     print_handler, 64932, (void *)(libmatlabTestInitializeWithHandlers)))
    return false;
  return true;
}

LIB_libmatlabTest_C_API 
bool MW_CALL_CONV libmatlabTestInitialize(void)
{
  return libmatlabTestInitializeWithHandlers(mclDefaultErrorHandler,
                                             mclDefaultPrintHandler);
}

LIB_libmatlabTest_C_API 
void MW_CALL_CONV libmatlabTestTerminate(void)
{
  if (_mcr_inst != NULL)
    mclTerminateInstance(&_mcr_inst);
}

LIB_libmatlabTest_C_API 
void MW_CALL_CONV libmatlabTestPrintStackTrace(void) 
{
  char** stackTrace;
  int stackDepth = mclGetStackTrace(_mcr_inst, &stackTrace);
  int i;
  for(i=0; i<stackDepth; i++)
  {
    mclWrite(2 /* stderr */, stackTrace[i], sizeof(char)*strlen(stackTrace[i]));
    mclWrite(2 /* stderr */, "\n", sizeof(char)*strlen("\n"));
  }
  mclFreeStackTrace(&stackTrace, stackDepth);
}


LIB_libmatlabTest_C_API 
bool MW_CALL_CONV mlxMatlabTest(int nlhs, mxArray *plhs[],
                                int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "matlabTest", nlhs, plhs, nrhs, prhs);
}

LIB_libmatlabTest_C_API 
bool MW_CALL_CONV mlfMatlabTest(int nargout, mxArray** vec_out, mxArray* vec_in)
{
  return mclMlfFeval(_mcr_inst, "matlabTest", nargout, 1, 1, vec_out, vec_in);
}

