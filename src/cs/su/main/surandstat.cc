/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SURANDSTAT: $Revision: 1.5 $ ; $Date: 2011/11/16 23:16:23 $	*/

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
#include "header.h"
extern "C" {
  #include <time.h>

}

/*********************** self documentation ******************************/
std::string sdoc_surandstat =
"									"
" SURANDSTAT - Add RANDom time shifts STATIC errors to seismic traces	"
"									"
"     surandstat <stdin >stdout  [optional parameters]	 		"
"									"
" Required parameters:							"
"	none								"
" Optional Parameters:							"
" 	seed=from_clock    	random number seed (integer)            "
"	max=tr.dt 		maximum random time shift (ms)		"
"	scale=1.0		scale factor for shifts			"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace surandstat {


/* Credits:
 *	U Houston: Chris Liner c. 2009
 */

/************************ end self doc ***********************************/


segy tr;

void* main_surandstat( void* args )
{
	int nt;		/* number of samples on output trace	*/
	float dt;	/* sample rate on outpu trace		*/
	int max;		/* max time shift (ms) */
	unsigned int seed;      /* random number seed */
	int it;			/* time sample counter */
	int itr;		/* trace count */
	int its;		/* local shift in time samples */
	float scale;		/* scale factor for shifts */


	/* Hook up getpar */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_surandstat );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get information from first trace */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nt   = tr.ns;
	dt   = ((double) tr.dt)/1000000.0;
	
	/* Get parameters */
	if (!parObj.getparint("max", &max))    max = (int) dt*1000.0;
	if (!parObj.getparfloat("scale", &scale))    scale = 1.0;


	/* Set seed */
	if (!parObj.getparuint("seed", &seed)) { /* if not supplied, use clock */
		if (-1 == (seed = (unsigned int) time((time_t *) NULL))) {
			throw cseis_geolib::csException("time() failed to set seed");
		}
	}
        parObj.checkpars();

	/* Loop on traces */	
	itr = 0;
	do {
		sranuni(seed);
		/* random shift (ms) for this trace */
		its = scale*2*max*(0.5 - franuni()) + 1;

		if (its <= 0 ) {
			/* loop over output times */
			for (it=0; it<nt-its; ++it) {
				tr.data[it] = tr.data[it-its];	
			}
		} else {
			/* loop over output times */
			for (it=nt; it>its; --it) {
				tr.data[it] = tr.data[it-its];	
			}

		}
		su2cs->putTrace(&tr);

		++itr;

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
