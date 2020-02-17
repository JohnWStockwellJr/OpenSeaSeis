/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUZERO: $Revision: 1.14 $ ; $Date: 2011/11/16 17:23:05 $	*/

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
std::string sdoc_suzero =
" 								"
" SUZERO -- zero-out (or set constant) data within a time window	"
" 								"
" suzero itmax= < indata > outdata				"
" 								"
" Required parameters						"
" 	itmax=		last time sample to zero out		"
" 								"
" Optional parameters						"
" 	itmin=0		first time sample to zero out		"
" 	value=0		value to set				"
" 								"
" See also: sukill, sumute					"
" 								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suzero {


/* Credits:
 *	CWP: Chris
 *	Geocon: Garry Perratt (added value= option)
 *
 * Trace header fields accessed: ns
 */
/**************** end self doc ***********************************/


segy tr;

void* main_suzero( void* args )
{
	int itmin;		/* first sample to zero out		*/
	int itmax;		/* last sample to zero out	 	*/
	float value;		/* value to set within window		*/
	int nt;			/* time samples per trace in input data	*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suzero );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get information from first trace */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;

	/* Get params from user */
	CSMUSTGETPARINT("itmax", &itmax);
	if (!parObj.getparint("itmin", &itmin))	itmin = 0;
	if (!parObj.getparfloat("value", &value))	value = 0.0;
        parObj.checkpars();

	/* Error checking */
	if (itmax > nt)    throw cseis_geolib::csException("itmax = %d, must be < nt", itmax);
	if (itmin < 0)     throw cseis_geolib::csException("itmin = %d, must not be negative", itmin);
	if (itmax < itmin) throw cseis_geolib::csException("itmax < itmin, not allowed");

	/* Main loop over traces */
	do { 
		register int i;
		for (i = itmin; i <= itmax; ++i)  tr.data[i] = value;
		
		su2cs->putTrace(&tr);
	} while(cs2su->getTrace(&tr));


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

} // END namespace
