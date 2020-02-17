/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUTSQ: $Revision: 1.8 $ ; $Date: 2011/11/16 23:21:55 $		*/

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
std::string sdoc_sutsq =
"								"
" SUTSQ -- time axis time-squared stretch of seismic traces	"
"								"
" sutsq [optional parameters] <stdin >stdout 			"
"								"
" Required parameters:						"
"	none				 			"
"								"
" Optional parameters:						"
"       tmin= .1*nt*dt  minimum time sample of interest		"
"                       (only needed for forward transform)	"
"       dt= .004       output sample rate			"
"                       (only needed for inverse transform)	"
"       flag= 1        1=forward transform: time to time squared"
"                     -1=inverse transform: time squared to time"
"								"
" Note: The output of the forward transform always starts with	"
" time squared equal to zero.  'tmin' is used to avoid aliasing	"
" the early times."
" "	
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sutsq {


/*
 * Caveats:
 * 	Amplitudes are not well preserved.
 *
 * Trace header fields accessed: ns, dt
 * Trace header fields modified: ns, dt
 */
/**************** end self doc ***********************************/

segy tr;

void* main_sutsq( void* args )
{
	int j,nt,flag,ntout;
	float *buf,*ttn,dt,dtout=0.0,tmin,tmax;

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sutsq );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get information from the first header */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;
	dt = (float) tr.dt/1000000.0;

	if (!parObj.getparfloat("tmin", &tmin)) tmin=0.1*nt*dt;
	if (!parObj.getparint("flag", &flag)) flag=1;
	if(flag==1) {
		dtout=tmin*2.*dt;
		tmax=nt*dt;
		ntout=1+tmax*tmax/dtout; CHECK_NT("ntout",ntout);
		ttn=ealloc1float(ntout);
		for(j=0;j<ntout;j++) ttn[j]=sqrt(j*dtout);
	}else{
		if (!parObj.getparfloat("dt", &dt)) dtout=0.004;
		ntout=1+sqrt(nt*dt)/dtout; CHECK_NT("ntout",ntout);
		ttn=ealloc1float(ntout);
		for(j=0;j<ntout;j++) ttn[j]=j*j*dtout*dtout;
	}
        parObj.checkpars();
	buf = ealloc1float(nt);

	fprintf(stderr,"sutsq: ntin=%d dtin=%f ntout=%d dtout=%f\n",
		nt,dt,ntout,dtout);

	/* Main loop over traces */
	do {
		for(j=0;j<nt;j++) buf[j]=tr.data[j];
		tr.ns = ntout;
		tr.dt = dtout*1000000.;			
		ints8r(nt,dt,0.,buf,0.0,0.0,
			ntout,ttn,tr.data);
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
