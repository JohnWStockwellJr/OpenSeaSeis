/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */


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
std::string sdoc_sucliphead =
" 									"
" SUCLIPHEAD - Clip header values					"
" 									"
" sucliphead <stdin >stdout [optional parameters]			"
"									"
" Required parameters:							"
"	none								"
"									"
" Optional parameters:							"
"	key=cdp,...			header key word(s) to clip	"
"	min=0,...			minimum value to clip		"
"	max=ULONG_MAX,ULONG_MAX,...	maximum value to clip		"
"									"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sucliphead {


/* Credits:
 *	Geocon: Garry Perratt
 *
 */
/**************** end self doc ********************************/

void changeval(cwp_String type, Value *valp, float fval);

segy tr;

void* main_sucliphead( void* args )
{
	cwp_String key[SU_NKEYS];	/* array of keywords		*/
	cwp_String type[SU_NKEYS];	/* array of types for key	*/
	int index[SU_NKEYS];		/* array of indexes for key	*/
	int ikey;		/* key counter				*/
	int nkeys;		/* number of header fields set		*/
	int n;			/* number of min,max values   		*/
	Value val;		/* value of key field			*/
	double fval;		/* value of key field			*/
	float *min=NULL;	/* array of "min" values		*/
	float *max=NULL;	/* array of "max" values		*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sucliphead );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get "key" values */
	if ((nkeys=parObj.countparval("key"))!=0) {
		parObj.getparstringarray("key",key);

	} else {
		key[0]="cdp";
	}

	/* get types and indexes corresponding to the keys */
	for (ikey=0; ikey<nkeys; ++ikey) {
		type[ikey]=hdtype(key[ikey]);
		index[ikey]=getindex(key[ikey]);
	}

	/* get "min" values */
	if ((n=parObj.countparval("min"))!=0) { 
		if (n!=nkeys)
		throw cseis_geolib::csException("number of a values not equal to number of keys");
		min=ealloc1float(n);
		parObj.getparfloat("min",min);
	} else {
		min=ealloc1float(nkeys);
		for (ikey=0; ikey<nkeys; ++ikey) min[ikey]=0.;
	}

	/* get "max" values */
	if ((n=parObj.countparval("max"))!=0) { 
		if (n!=nkeys)
		throw cseis_geolib::csException("number of a values not equal to number of keys");
		max=ealloc1float(n);
		parObj.getparfloat("max",max);
	} else {
		max=ealloc1float(nkeys);
		for (ikey=0; ikey<nkeys; ++ikey) max[ikey]=ULONG_MAX;
	}

        parObj.checkpars();
	/* get types and index values */
	for (ikey=0; ikey<nkeys; ++ikey) {
		type[ikey] = hdtype(key[ikey]);
		index[ikey] = getindex(key[ikey]);
	}

	while (cs2su->getTrace(&tr)) {
		for (ikey=0; ikey<nkeys; ++ikey) {
			gethval(&tr, index[ikey], &val);
			fval = vtof(type[ikey], val);
			if (fval < min[ikey]) {
				changeval(type[ikey], &val, min[ikey]);
				puthval(&tr, index[ikey], &val);
			} else if (fval > max[ikey]) {
				changeval(type[ikey], &val, max[ikey]);
				puthval(&tr, index[ikey], &val);
			}
		}
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

void changeval(cwp_String type, Value *valp, float fval)
{
	switch (*type) {
	case 's':
		throw cseis_geolib::csException("cannot change char header word");
	break;
	case 'h':
		valp->h = (short) fval;
	break;
	case 'u':
		valp->u = (unsigned short) fval;
	break;
	case 'l':
		valp->l = (long) fval;
	break;
	case 'v':
		valp->v = (unsigned long) fval;
	break;
	case 'i':
		valp->i = (int) fval;
	break;
	case 'p':
		valp->p = (unsigned int) fval;
	break;
	case 'f':
		valp->f = (float) fval;
	break;
	case 'd':
		valp->d = (double) fval;
	break;
	default:
		throw cseis_geolib::csException("unknown type %s", type);
	break;
	}
}

} // END namespace
