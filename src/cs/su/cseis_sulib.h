/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CSEIS_SU_LIB_H
#define CSEIS_SU_LIB_H

#include "csSUGetPars.h"

/* Forward declarations of SU library functions which had to be modified for SEASEIS  */

void cseis_decodeReflectors( cseis_su::csSUGetPars* parObj, int *nrPtr,
                             float **aPtr, int **nxzPtr, float ***xPtr, float ***zPtr);

void cseis_getParVal( cseis_su::csSUGetPars* parObj, cwp_String name, cwp_String type, int n, Value *valp);

#endif
