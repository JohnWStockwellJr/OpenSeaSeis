/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUDIVCOR: $Revision: 1.15 $ ; $Date: 2011/11/16 17:23:05 $		*/

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

/*********************** self documentation ******************************/
std::string sdoc_sudivcor =
" 									"
" SUDIVCOR - Divergence (spreading) correction				"
" 									"
" sudivcor <stdin >stdout  [optional parms]				"
" 									"
" Required Parameters:							"
" none									"
" 									"
" Optional Parameters:							"
" trms=0.0	times corresponding to rms velocities in vrms		"
" vrms=1500.0	interval velocities corresponding to times in trms	"
" vfile=	binary (non-ascii) file containing velocities vrms(t)	"
" 									"
" Notes:								"
" The trms, vrms arrays specify an rms velocity function of time.	"
" Linear interpolation and constant extrapolation is used to determine	"
" interval velocities at times not specified.  Values specified in trms	"
" must increase monotonically.						"
" 									"
" Alternatively, rms velocities may be stored in a binary file		"
" containing one velocity for every time sample.  If vfile is specified,"
" then the trms and vrms arrays are ignored.				"
" 									"
" The time of the first sample is assumed to be constant, and is taken	"
" as the value of the first trace header field delrt. 			"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sudivcor {


/* Credits:
 *	CWP: Jack K. Cohen, Francesca Fazarri
 *
 * Trace header fields accessed:  ns, dt, delrt
 */
/**************** end self doc *******************************************/


segy tr;

void* main_sudivcor( void* args )
{
	int nt;			/* number of points on input trace	*/
	float dt;		/* sample spacing			*/
	float tmin;		/* time of first sample			*/
	float *vt;		/* velocity function			*/
	float *vrms;		/* rms velocity picks			*/
	float *trms;		/* times corresponding to vrms picks	*/
	cwp_String vfile="";	/* binary file giving vrms(t)		*/
	float *divcor;		/* divergence correction function	*/
	float denom=1.0; 	/* vrms to divide by 			*/

	/* Hook up getpar */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sudivcor );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get information from the first header */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;
	dt = ((double) tr.dt)/1000000.0;
	tmin = tr.delrt/1000.0;


	/* Determine velocity function v(t) */
	vt = ealloc1float(nt);
	if (!parObj.getparstring("vfile", &vfile)) {
		int ntrms = parObj.countparval("trms");
		int nvrms = parObj.countparval("vrms");
		int npar;
		register int it, itrms;

		if (nvrms != ntrms)
			throw cseis_geolib::csException("number of trms and vrms must be equal");

		if (!ntrms)  ntrms = 2;  /* default case */
		trms = ealloc1float(ntrms);
		if (!(npar = parObj.getparfloat("trms", trms))) {
			trms[0] = 0.0;
			trms[1] = dt;
		}
		if (npar == 1) trms[1] = trms[0] + dt; /* const vel case */
		
		if (!nvrms)  nvrms = 2;  /* default case */
		vrms = ealloc1float(nvrms);
		if (!(npar = parObj.getparfloat("vrms", vrms))) {
			vrms[0] = 1500.0;
			vrms[1] = 1500.0;
		}
		if (npar == 1) vrms[1] = vrms[0];  /* const vel case */

		for (itrms = 1; itrms < ntrms; ++itrms)
			if (trms[itrms] <= trms[itrms-1])
				throw cseis_geolib::csException("trms must increase monotonically");

		for (it = 0; it < nt; ++it) {
			float t = tmin + it*dt;
			intlin(ntrms,trms,vrms,vrms[0],vrms[ntrms-1],
				1,&t,&vt[it]);
		}
		denom = trms[1]*vrms[1]*vrms[1];
		
	} else {  /* user gave a vfile */
		FILE *fp = efopen(vfile, "r");
		
		if (nt != efread(vt, FSIZE, nt, fp)) {
			throw cseis_geolib::csException("cannot read %d velocities from file %s",
				nt, vfile);
		} else denom = (tmin + dt)*vt[1]*vt[1];
	}
        parObj.checkpars();

	/* Form divergence correction vector */
	{ register int it;
	  
	  divcor = ealloc1float(nt);
	  for (it = 0; it < nt; ++it) {
	  	float t = tmin + it*dt;
		divcor[it] = t*vt[it]*vt[it] / denom;
	  }
	}
   
	  
	/* Main loop over traces */
	do {
		register int it;

	  	for (it = 0; it < nt; ++it)  tr.data[it] *=  divcor[it];
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
