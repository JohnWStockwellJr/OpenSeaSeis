/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUABSHW: $Revision: 1.13 $ ; $Date: 2011/11/16 22:10:29 $	*/

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
std::string sdoc_suabshw =
" 								"
" SUABSHW - replace header key word by its absolute value	"
" 								"
" suabshw <stdin >stdout key=offset				"
" 								"
" Required parameters:						"
" 	none							"
" 								"
" Optional parameter:						"
" 	key=offset		header key word			"
" 								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suabshw {


/* Credits:
 *	CWP: Jack K. Cohen
 */
/**************** end self doc ***********************************/


#define KEY	"offset"	/* Default key word to take abs() of */

/* function prototype of function used internally */
void absval(cwp_String type, Value *valp);

segy tr;

void* main_suabshw( void* args )
{
	cwp_String key;	/* header key word from segy.h		*/
	cwp_String type;/* ... its type				*/
	int index;	/* ... its index in hdr.h		*/
	Value val;	/* ... its value			*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suabshw );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get key parameter */
	if (!parObj.getparstring("key", &key))	key = KEY;

        parObj.checkpars();

	type = hdtype(key);
	index = getindex(key);

	while (cs2su->getTrace(&tr)) {
		gethval(&tr, index, &val);
		absval(type, &val);
		puthval(&tr, index, &val);
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


void absval(cwp_String type, Value *valp)
{
	switch (*type) {
	case 's': throw cseis_geolib::csException("can't absval char header word"); break;
	case 'u':	/* do nothing if unsigned type */
	case 'v':
	case 'p':                                      break;
	case 'h': if (valp->h < 0) valp->h = -valp->h; break;
	case 'l': if (valp->l < 0) valp->l = -valp->l; break;
	case 'i': if (valp->i < 0) valp->i = -valp->i; break;
	case 'f': if (valp->f < 0) valp->f = -valp->f; break;
	case 'd': if (valp->d < 0) valp->d = -valp->d; break;
	default: throw cseis_geolib::csException("unknown type %s", type);         break;
	}
}

} // END namespace
