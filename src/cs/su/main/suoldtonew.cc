/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUOLDTONEW: $Revision: 1.5 $ ; $Date: 2011/11/12 00:01:45 $	*/

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

/*********************** self documentation **********************/
std::string sdoc_suoldtonew =
" 								"
" SUOLDTONEW - convert existing su data to xdr format		"
" 								"
" suoldtonew <oldsu >newsu  					"
"								"
" Required parameters:						"
"	none							"
" 								"
" Optional parameters:						"
"	none							"
"								"
" Notes:							"
" This program is used to convert native machine datasets to	"
" xdr-based, system-independent format.				"
"								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suoldtonew {


/*
 * Author: Stewart A. Levin, Mobil, 1966
 *  
 */
/**************** end self doc ***********************************/


segy tr;

void* main_suoldtonew( void* args )
{

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suoldtonew );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	while (!(feof(stdin) || ferror(stdin))) {
		static int ntr=0; /* for user info only */

		/* Do read of header for the segy */
		if (0 >= efread(&tr, HDRBYTES, 1, stdin)) {
			warn("converted %d traces to XDR format", ntr);
			break; /* loop exit */
		}

		/* Do read of data for the segy */
		switch(efread((char *) (&(tr.data[0])), FSIZE, tr.ns, stdin)) {
		case 0: /* oops, no data for this header */
			throw cseis_geolib::csException("header without data for trace #%d", ntr+1);
		default:
			su2cs->putTrace(&tr);
			++ntr;
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
