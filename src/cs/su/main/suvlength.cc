/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUVLENGTH: $Revision: 1.8 $ ; $Date: 2011/11/16 23:09:52 $	*/

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
std::string sdoc_suvlength =
" 								"
" SUVLENGTH - Adjust variable length traces to common length   	"
" 								"
" suvlength <vdata >stdout					"
" 								"
" Required parameters:						"
" 	none							"
" 								"
" Optional parameters:						"
" 	ns=	output number of samples (default: 1st trace ns)"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suvlength {

/**************** end self doc ***********************************/

/* Credits:
 *	CWP: Jack Cohen, John Stockwell
 *
 * Trace header fields accessed:  ns
 * Trace header fields modified:  ns
 */

segy tr;

void* main_suvlength( void* args )
{
	int ns;		/* samples on output traces	*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suvlength );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



        /* Get info from first trace */ 
        if (!fvcs2su->getTrace(stdin, &tr))  throw cseis_geolib::csException("can't get first trace"); 
        if (!parObj.getparint("ns", &ns)) ns = tr.ns;
        parObj.checkpars();


	/* Loop over the traces */
	do {
		int nt = tr.ns;
				
		if (nt < ns) /* pad with zeros */
                	memset((void *)(tr.data + nt), 0, (ns-nt)*FSIZE);
		tr.ns = ns;
		su2cs->putTrace(&tr);
	} while (fvcs2su->getTrace(stdin, &tr));
	
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
