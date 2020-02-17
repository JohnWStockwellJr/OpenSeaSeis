/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUREDUCE: $Revision: 1.13 $ ; $Date: 2011/11/16 23:21:55 $	*/

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
std::string sdoc_sureduce =
" 	   							"
" SUREDUCE - convert traces to display in reduced time		" 
" 	   							"
" sureduce <stdin >stdout rv=					"
" 								"
" Required parameters:						"
"	dt=tr.dt	if not set in header, dt is mandatory	"
"								"
" Optional parameters:						"
"	rv=8.0		reducing velocity in km/sec		"	
"								"
" Note: Useful for plotting refraction seismic data. 		"
" To remove reduction, do:					"
" suflip < reduceddata.su flip=3 | sureduce rv=RV > flip.su	"
" suflip < flip.su flip=3 > unreduceddata.su			"
"								"
" Trace header fields accessed: dt, ns, offset			"
" Trace header fields modified: none				"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sureduce {


/*
 * Author: UC Davis: Mike Begnaud  March 1995
 *
 *
 * Trace header fields accessed: ns, dt, offset
 */
/**************** end self doc ***********************************/
segy tr;

void* main_sureduce( void* args )
{
	int nt;
	float dt, rv;

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sureduce );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get parameters */
	if (!parObj.getparfloat("rv", &rv))	 rv=8.0;

	/* Data validation */
	if (rv <= 0.0) throw cseis_geolib::csException("rv=%f, must be >= 0.0", rv);

	/* Get info from first trace */
	if (!cs2su->getTrace( &tr)) throw cseis_geolib::csException("can't read first trace");
	if (!tr.dt) CSMUSTGETPARFLOAT("dt", &dt);
	nt = (int) tr.ns;
	dt = ((double) tr.dt)/1000000.0;	
        parObj.checkpars();

	/* Loop over traces */
	do {
		float off  = (float) ABS(tr.offset);
		float bt   = ((fabs) (off))/(rv*1000.0);
		int rnt    = NINT(bt/dt);
		register int i;

		for (i = 0; i < nt; ++i) {
			register int j = i + rnt;
			tr.data[i] = (i < (nt - rnt)) ?  tr.data[j] : 0.0;
		}
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


} // END namespace
