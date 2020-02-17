/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SEGYCLEAN: $Revision: 1.13 $ ; $Date: 2011/11/12 00:01:45 $		*/

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
std::string sdoc_segyclean =
" 								"
" SEGYCLEAN - zero out unassigned portion of header		"
" 								"
" segyclean <stdin >stdout 					"
"								"
" Since \"foreign\" SEG-Y tapes may use the unassigned portion	"
" of the trace headers and since SU now uses it too, this	"
" program zeros out the fields meaningful to SU.		"
" 								"
"  Example:							"
"  	segyread trmax=200 | segyclean | suximage		"
"								"
"								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace segyclean {


/* Credits:
 *	CWP: Jack Cohen
 *
 */
/**************** end self doc ********************************/


segy tr;

void* main_segyclean( void* args )
{

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_segyclean );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	while (cs2su->getTrace(&tr)) {
		tr.f1 = 0.0;
		tr.d1 = 0.0;
		tr.f2 = 0.0;
		tr.d2 = 0.0;
		tr.ungpow = 0.0;
		tr.unscale = 0.0;
		tr.ntr = 0;
		tr.mark = 0;

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
