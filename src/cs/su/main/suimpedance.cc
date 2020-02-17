/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUIMPEDANCE: $Revision: 1.4 $ ; $Date: 2011/11/12 00:40:42 $	*/

#include "csException.h"
#include "csSUTraceManager.h"
#include "csSUArguments.h"
#include "csSUGetPars.h"
#include "su_complex_declarations.h"
#include "cseis_sulib.h"
#include <string>

extern "C" {
  #include <pthread.h>
}
#include "su.h"
#include "segy.h"


/*********************** self documentation **********************/
std::string sdoc_suimpedance =
"								"
" SUIMPEDANCE - Convert reflection coefficients to impedances.  "
"								"
" suimpedance <stdin >stdout [optional parameters]		"
"								"
" Optional Parameters:					  	"
" v0=1500.	Velocity at first sample (m/sec)		"
" rho0=1.0e6	Density at first sample  (g/m^3)		"
"								"
" Notes:							"
" Implements recursion [1-R(k)]Z(k) = [1+R(k)]Z(k-1).		"
" The input traces are assumed to be reflectivities, and thus are"
" expected to have amplitude values between -1.0 and 1.0.	"
"								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suimpedance {


/* Credits:
 *	SEP: Stew Levin
 *
 * Trace header fields accessed: ns
 * 
 */
/**************** end self doc ***********************************/

static segy tr;

/* Prototype of function used internally */
static void rctoimp(float v0, float rho0, int nt, float *trace);

void* main_suimpedance( void* args )
{
  float v0, rho0;
  int nt;
  
  /* Initialize */
  cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
  cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
  cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
  int argc = suArgs->argc;
  char **argv = suArgs->argv;
  cseis_su::csSUGetPars parObj;

  void* retPtr = NULL;  /*  Dummy pointer for return statement  */
  su2cs->setSUDoc( sdoc_suimpedance );
  if( su2cs->isDocRequestOnly() ) return retPtr;
  parObj.initargs(argc, argv);

  try {  /* Try-catch block encompassing the main function body */

  
  /* set parameters and fill header fields */
  if (!parObj.getparfloat("v0", &v0)) v0 = 1500.;
  if (!parObj.getparfloat("rho0", &rho0)) rho0 = 1000000.;

  /* Get info from first trace */
  if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
  nt = tr.ns;

  do {
	/* apply reflectivity to impedence operation */
	rctoimp( v0, rho0, nt, tr.data );
	su2cs->putTrace(&tr); 	 	

	} while (cs2su->getTrace(&tr));

	su2cs->setEOF();
	pthread_exit(NULL);
	return retPtr;
}
catch( cseis_geolib::csException& exc ) {
  su2cs->setError("%s",exc.getMessage());
  pthread_exit(NULL);
  return retPtr;
}
}

static void
rctoimp(float v0, float rho0, int nt, float *trace)
/***************************************************************************
rctoimp - convert reflection coefficient trace to impedence trace
****************************************************************************
Author: SEP: Stew Levin
**************************************************** ***********************/
{
  int it;
  double  zold, rc;

  zold = v0*rho0;

  for(it = 0; it < nt; ++it) {
	rc = trace[it];
	if(rc <= -1.0) rc = -0.9999;
	if(rc >=  1.0) rc =  0.9999;
	trace[it] = (float) zold;
	zold *= (1.0+rc) / (1.0-rc);
  }
}

} // END namespace
