/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUNULL: $Revision: 1.15 $ ; $Date: 2011/11/12 00:40:42 $	*/

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
std::string sdoc_sunull =
" 								"
" SUNULL - create null (all zeroes) traces	 		"
" 								"
" sunull nt=   [optional parameters] >outdata			"
" 								"
" Required parameter						"
" 	nt=		number of samples per trace		"
" 								"
" Optional parameters						"
" 	ntr=5		number of null traces to create		"
" 	dt=0.004	time sampling interval			"
" 								"
" Rationale: It is sometimes useful to insert null traces	"
"	 between \"panels\" in a shell loop.			"
" 								"
" See also: sukill, sumute, suzero				"
" 								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sunull {


/* Credits:
 *	CWP: Jack K. Cohen
 *
 * Trace header fields set: ns, dt, tracl
 */
/**************** end self doc ***********************************/


segy tr;

void* main_sunull( void* args )
{
	int nt;			/* number of time samples		*/
	int ntr;		/* number of traces			*/
	register int itr;	/* trace counter			*/
	float dt;		/* time sampling interval (seconds)	*/
	int idt;		/* 	...		(micro seconds)	*/
	
	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sunull );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get parameters */
	CSMUSTGETPARINT("nt", &nt);
	CHECK_NT("nt",nt);
        tr.ns = nt;

	if (!parObj.getparint("ntr", &ntr))	ntr = 5;
	if (!parObj.getparfloat("dt", &dt))	dt = .004;
	idt = 1000000.0 * dt;

	/* Set tr.data to zeros */
	memset( (void *) tr.data, 0, nt*FSIZE);

	/* Set constant header fields */
	tr.dt = idt;
	tr.ns = nt;

	/* Main loop over traces */
	for (itr = 0; itr < ntr; ++itr) {

		tr.tracl = itr + 1;

		su2cs->putTrace(&tr);
	}

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
