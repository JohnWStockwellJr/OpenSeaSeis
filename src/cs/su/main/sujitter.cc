/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUJITTER: $Revision: 1.3 $ ; $Date: 2011/11/12 00:22:43 $	*/

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
std::string sdoc_sujitter =
"									"
" SUJITTER - Add random time shifts to seismic traces			"
"									"
"     sujitter <stdin >stdout  [optional parameters]	 		"
"									"
" Required parameters:							"
"	none								"
" Optional Parameters:							"
" 	seed=from_clock    	random number seed (integer)            "
"	min=1 			minimum random time shift (samples)	"
"	max=1 			maximum random time shift (samples)	"
"	pon=1 			shift can be positive or negative	"
"				=0 shift is positive only		"
"	fldr=0 			each trace has new shift		"
"				=1 new shift when fldr header field changes"
" Notes:								"
" Useful for simulating random statics. See also:  suaddstatics		"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sujitter {


/* Credits:
 *	U of Houston: Chris Liner 
 *	UH:  Chris added fldr, min, pon options 12/10/08
 */

/************************ end self doc ***********************************/


segy tr;

void* main_sujitter( void* args )
{
	int nt;		/* number of samples on output trace	*/
	float dt;	/* sample rate on outpu trace		*/
	int min;		/* min time shift (samples) */
	int max;		/* max time shift (samples) */
	unsigned int seed;      /* random number seed */
	int it;			/* time sample counter */
	int itr;		/* trace counter */
	int its=0;		/* local shift in time samples */
	int sits;		/* sign of local shift */
	int fldr;		/* fldr use flag		*/
	int ishot;		/* shot counter (based on tr.fldr) */
	int ishotold;
	int pon;		/* flag for pos or neg shift */


	/* Hook up getpar */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sujitter );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get information from first trace */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nt   = tr.ns;
	ishotold = tr.fldr;
	dt   = ((double) tr.dt)/1000000.0;
	
	/* Get parameters */
	if (!parObj.getparint("min", &min))    min = 1;
	if (!parObj.getparint("max", &max))    max = 1;
	if (!parObj.getparint("pon", &pon))    pon = 1;
	if (!parObj.getparint("fldr", &fldr))  fldr = 0;

	if (min>max) throw cseis_geolib::csException("min>max... exit");

	/* Set seed */
	if (!parObj.getparuint("seed", &seed)) { /* if not supplied, use clock */
		if (-1 == (seed = (unsigned int) time((time_t *) NULL))) {
			throw cseis_geolib::csException("time() failed to set seed");
		}
	}
	sranuni(seed);

	/* Loop on traces */	
	itr = 1;
	do {

		/* get fldr (shot) number */
		ishot = tr.fldr;

		if (itr==1) {
			/* initial shift for shot 1 
			   (used if fldr==1) */	
	       		its = min + (max-min)*franuni();
			if (pon==1) {
				/* include random sign to shift */
				sits = SGN(franuni()-0.5);
				its = its * sits;
			}
		}
		
		if (fldr==0) {
			/* each trace gets random shift */
	       		its = min + (max-min)*franuni();
			if (pon==1) {
				/* include random sign to shift */
				sits = SGN(franuni()-0.5);
				its = its * sits;
			}
		}

		if (fldr==1 && ishot!=ishotold) {
			/* new shot needs new shift */
	       		its = min + (max-min)*franuni();
			if (pon==1) {
				/* include random sign to shift */
				sits = SGN(franuni()-0.5);
				its = its * sits;
			}
			ishotold = ishot;
		}

		/* apply shift and output trace */
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

		itr += 1;

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
