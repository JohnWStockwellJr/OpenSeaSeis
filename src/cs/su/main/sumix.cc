/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUMIX: $Revision: 1.11 $ ; $Date: 2011/11/16 23:09:52 $	*/

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
std::string sdoc_sumix =
"									"
" SUMIX - compute weighted moving average (trace MIX) on a panel	"
"	  of seismic data						"
"									"
" sumix <stdin >sdout 							"
" mix=.6,1,1,1,.6	array of weights for weighted average		"
"									"
"									"
" Note: 								"
" The number of values defined by mix=val1,val2,... determines the number"
" of traces to be averaged, the values determine the weights.		"
" 									"
" Examples: 								"
" sumix <stdin mix=.6,1,1,1,.6 >sdout 	(default) mix over 5 traces weights"
" sumix <stdin mix=1,1,1 >sdout 	simple 3 trace moving average	"
" 									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sumix {


/* Author:
 *	CWP: John Stockwell, Oct 1995
 *
 * Trace header fields accessed: ns
 */
/**************** end self doc ***********************************/


/* default weighting values */
#define VAL0	0.6
#define VAL1	1.0
#define VAL2	1.0
#define VAL3	1.0
#define VAL4	0.6


segy tr;

void* main_sumix( void* args )
{
	int nmix;		/* number of traces to mix over		*/
	int imix;		/* mixing counter			*/
	int it;			/* sample counter			*/
	int nt;			/* number of time samples per trace	*/
	int itr=0;		/* trace counter			*/
	size_t databytes;	/* number of bytes (nt*FSIZE)		*/
	size_t mixbytes;	/* number of bytes (nt*FSIZE*nmix)	*/
	float *mix;		/* array of mix values			*/
	float *temp;		/* temp array for mixing 		*/
	float **data;		/* array for mixing 			*/
	
	
	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sumix );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get info from first trace */
	if(!cs2su->getTrace(&tr))
		throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;

	/* Get mix weighting values values */
	if ((nmix = parObj.countparval("mix"))!=0) {
		mix = ealloc1float(nmix);
		parObj.getparfloat("mix",mix);
		
	} else {
		nmix = 5;
		mix = ealloc1float(nmix);

		mix[0] = VAL0;
		mix[1] = VAL1;
		mix[2] = VAL2;
		mix[3] = VAL3;
		mix[4] = VAL4;
	}

        parObj.checkpars();

	/* Divide mixing weight by number of traces to mix */
	for (imix = 0; imix < nmix; ++imix)
		mix[imix]=mix[imix]/((float) nmix);

	/* Compute databytes per trace and bytes in mixing panel */
	databytes = FSIZE*nt;
	mixbytes = databytes*nmix;

	/* Allocate temporary space for mixing  */
	data = ealloc2float(nt,nmix);
	temp = ealloc1float(nt);

	/* Zero out data array */
	memset((void *) data[0], 0, mixbytes);

	/* Loop over remaining traces */
	do {

		++itr;

		/* Zero out temp */
		memset((void *) temp, 0, databytes);

		/* Read data portion of trace into first column of data[][] */
		memcpy( (void *) data[0], (const void *) tr.data, databytes);
	
		/* Loop over time samples */
		for (it=0; it<nt; ++it) {

			/* Weighted moving average (mix) */
			for(imix=0; imix<nmix; ++imix)
				temp[it]+=data[imix][it]*mix[imix];

			/* put mixed data back in seismic trace */
			tr.data[it] = temp[it]; 
		}

		/* Bump columns of data[][] over by 1 */
		/* to make space for data from next trace */
		for (imix=nmix-1; 0<imix; --imix)
			for (it=0; it<nt; ++it) 
				data[imix][it] = data[imix-1][it];

				
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
