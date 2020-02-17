/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUACOR: $Revision: 1.15 $ ; $Date: 2011/11/16 17:37:27 $		*/

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

/*********************** self documentation ******************************/
std::string sdoc_suacor =
"									"
" SUACOR - auto-correlation						"
"									"
" suacor <stdin >stdout [optional parms]				"
"									"
" Optional Parameters:							"
" ntout=101	odd number of time samples output			"
" norm=1	if non-zero, normalize maximum absolute output to 1	"
" sym=1		if non-zero, produce a symmetric output from		"
"			lag -(ntout-1)/2 to lag +(ntout-1)/2		"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suacor {


/* Credits:
 *	CWP: Dave Hale
 *
 * Trace header fields accessed:  ns
 * Trace header fields modified:  ns and delrt
 */
/**************** end self doc *******************************************/

segy tr;

void* main_suacor( void* args )
{
	int nt,ntout,it,istart,izero,norm,sym;
	float scale,*temp;

	/* hook up getpar */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suacor );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* get information from the first header */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;

	/* get parameters */
	if (!parObj.getparint("ntout",&ntout)) ntout=101;
	if (!parObj.getparint("norm",&norm)) norm = 1;
	if (!parObj.getparint("sym",&sym)) sym = 1;
        parObj.checkpars();
	
	/* allocate workspace */
	temp = ealloc1float(ntout);
	
	/* index of first sample */
	if (sym == 0) istart = 0;
	else istart = -(ntout-1)/2;

	/* index of sample at time zero */
	izero = -istart;
	
	/* loop over traces */
	do {
		xcor(nt,0,tr.data,nt,0,tr.data,ntout,istart,temp);
		if (norm) {
			scale = 1.0/(temp[izero]==0.0?1.0:temp[izero]);
			for (it=0; it<ntout; ++it)  temp[it] *= scale;
		}
		memcpy((void *) tr.data, (const void *) temp, ntout*FSIZE);
		tr.ns = ntout;
		tr.delrt = 0;
		su2cs->putTrace(&tr);
	} while(cs2su->getTrace(&tr));

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
