/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUHILB: $Revision: 1.19 $ ; $Date: 2011/11/12 00:42:19 $	*/

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
std::string sdoc_suhilb =
"								"
" SUHILB - Hilbert transform					"
"								"
" suhilb <stdin >sdout 						"
"							        "
" Note: the transform is computed in the direct (time) domain   "
"							        "
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suhilb {


/* Credits:
 *	CWP: Jack Cohen   
 *      CWP: John Stockwell, modified to use Dave Hale's hilbert() subroutine
 *
 * Trace header fields accessed: ns, trid
 */
/**************** end self doc ***********************************/


segy tr;

void* main_suhilb( void* args )
{
	int nt;			/* number of points on input trace	*/
	float *data;		/* data values from each trace		*/
	float *hdata;		/* Hilbert transformed data values	*/
	register int i;		/* counter				*/
	cwp_Bool seismic;	/* is this seismic data?		*/
	
	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suhilb );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get info from first trace */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;

	/* Check that data is correct format */
	seismic = ISSEISMIC(tr.trid);
	if (!seismic)
		warn("input is not seismic data, trid=%d", tr.trid);

	/* allocate space for data and hdata */
	data = ealloc1float(nt);
	hdata = ealloc1float(nt);

	/* Loop over traces */
	do {

		/* read data from trace */
		for (i = 0; i < nt; i++) data[i] = tr.data[i];

		/* apply the Hilbert tranform */		
		hilbert(nt,data,hdata);

		/* put Hilbert tranformed data back in trace */
		for (i = 0; i < nt; i++) tr.data[i] = hdata[i];
		
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
