/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SWAPBHED: $Revision: 1.2 $ ; $Date: 2011/11/12 00:01:45 $	*/

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
extern "C" {
  #include <stdint.h>

}
#include "su.h"
#include "segy.h"
#include "bheader.h"

/*********************** self documentation **********************/
std::string sdoc_swapbhed =
" 									"
" SWAPBHED - SWAP the BYTES in a SEGY Binary tape HEaDer file		"
" 									"
" swapbhed < binary_in > binary out					"
" 									"
" Required parameter:							"
" 	none								"
" Optional parameters:							"
"	none 								"
" 									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace swapbhed {


/* Credits:
 *
 *	CWP: John Stockwell  13 May 2011
 */
/**************** end self doc ***********************************/

segy tr;
bhed bh;

void* main_swapbhed( void* args )
{
	int i;			/* counter				*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_swapbhed );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Read in binary header from standard in  */
	efread(&bh,BNYBYTES, 1, stdin);

	/* swap bytes */
	for (i = 0; i < BHED_NKEYS; ++i) swapbhval(&bh, i);

	/* Write binary header from bh structure to standard out */
	efwrite( (char *) &bh, 1, BNYBYTES, stdout);

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
