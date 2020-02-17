/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUKILL: $Revision: 1.18 $ ; $Date: 2011/11/17 00:03:38 $	*/

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
std::string sdoc_sukill =
" 								"
" SUKILL - zero out traces					"
" 								"
" sukill <stdin >stdout [optional parameters]			"
" 								"
" Optional parameters:						"
"	key=trid	header name to select traces to kill	"
"	a=2		header value identifying tracces to kill"
" or								"
" 	min= 		first trace to kill (one-based)		"
" 	count=1		number of traces to kill 		"
" 								"
" Notes:							"
"	If min= is set it overrides selecting traces by header.	"
" 								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sukill {


/* Credits:
 *	CWP: Chris Liner, Jack K. Cohen
 *	header-based trace selection: Florian Bleibinhaus
 *
 * Trace header fields accessed: ns
 */
/**************** end self doc ***********************************/


segy tr;

void* main_sukill( void* args )
{
	cwp_String key;		/* trace header			*/
	cwp_String type;	/* type for trace header	*/
	int index;		/* index of trace header	*/
	Value val;		/* trace header value		*/
	double dval,a;		/* trace header value		*/
	register int itr;	/* trace counter		*/
	int min;		/* first trace to zero out	*/
	int count;		/* number of traces to zero out	*/
	int nt = 0;		/* number of time samples	*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sukill );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get parameters */
	if (!parObj.getparint("min", &min)) min = 0;
	if (min < 0) throw cseis_geolib::csException("min = %d, must be >= 1", min);
	if (!parObj.getparint("count", &count))  count = 1;
	if (!parObj.getparstring("key", &key)) key = "trid";
	if (!parObj.getpardouble("a", &a)) a = 2.;
        parObj.checkpars();

	/* Get type and index value */
	type  = hdtype(key);
	index = getindex(key);


	if ( min>0 ) {

	/* Echo traces till min if they are there */
	for (itr = 1; itr < min; ++itr) {
		if (cs2su->getTrace(&tr)) su2cs->putTrace(&tr);
		else throw cseis_geolib::csException("failed to get requested trace #%ld", itr);
	}

	/* Kill "count" traces if they are there
	 * Do first outside loop to get nt    */
	if (cs2su->getTrace(&tr)) {
		nt = tr.ns;
		memset( (void *) tr.data, 0, nt * FSIZE);
		su2cs->putTrace(&tr);
		++itr;
	} else throw cseis_geolib::csException("failed to get requested trace #%ld", itr);

	for ( ; itr < min + count; ++itr) {
		if (cs2su->getTrace(&tr)) {
			memset( (void *) tr.data, 0, nt * FSIZE);
			su2cs->putTrace(&tr);
		} else throw cseis_geolib::csException("failed to get requested trace #%ld", itr);
	}

	/* Echo the trailing traces if any */
	while (cs2su->getTrace(&tr)) {
		su2cs->putTrace(&tr);
	}

	} else {	/* select traces by header value */
		while (cs2su->getTrace(&tr)) {
			nt = tr.ns;
			gethval(&tr, index, &val);
			dval = vtod(type, val);
			if ( dval==a ) memset( (void *) tr.data, 0, nt*FSIZE);
			su2cs->putTrace(&tr);
		}
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
