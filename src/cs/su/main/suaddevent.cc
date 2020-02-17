/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.		       */

/* SUADDEVENT: $Revision: 1.6 $ ; $Date: 2011/11/16 23:30:27 $	       */

#define SQ(x) ((x))*((x))
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

std::string sdoc_suaddevent =
"								       "
" SUADDEVENT - add a linear or hyperbolic moveout event to seismic data "
"								       "
" suaddevent <stdin >stdout [optional parameters]		       "
"								       "
" Required parameters:						  "
"       none								"
"								       "
" Optional parameters:						  "
"     type=nmo    =lmo for linear event 				"
"     t0=1.0      zero-offset intercept time IN SECONDS			"
"     vel=3000.   moveout velocity in m/s				"
"     amp=1.      amplitude						"
"     dt=	 must provide if 0 in headers (seconds)		"
"   									"
" Typical usage: "
"     sunull nt=500 dt=0.004 ntr=100 | sushw key=offset a=-1000 b=20 \\ "
"     | suaddevent v=1000 t0=0.05 type=lmo | suaddevent v=1800 t0=0.8 \\"
"     | sufilter f=8,12,75,90 | suxwigb clip=1 &	     		"
"								       "
"								       "
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suaddevent {


/* Credits:
 *      Gary Billings, Talisman Energy, May 1996, Apr 2000, June 2001
 *
 * Note:  code is inefficient in that to add a single "spike", with sinc
 *	interpolation, an entire trace is generated and added to 
 *	the input trace.  In fact, only a few points needed be created
 *	and added, but the current coding avoids the bookkeeping re
 *	which are the relevant points!
 */
/**************** end self doc *******************************************/

segy tr;

void* main_suaddevent( void* args )
{

	int nsegy, nt;
	int i;
	float dt, tmin, t0, vel, tx;
	float amp;
	char *type;
	int typecode=0;
	float *addtr;
	float *times;

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suaddevent );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	
	nsegy=cs2su->getTrace(&tr);
	if (!nsegy) throw cseis_geolib::csException("cannot get first trace");

	/* Get nt, etc. from first trace */
	nt   = (int) tr.ns;
	dt = ((double) tr.dt)/1000000.0;   /* microsecs to secs */
	tmin = tr.delrt/1000.0;		   /* millisecs to secs */

	addtr=alloc1float( nt );
	times=alloc1float( nt );
	for(i=0;i<nt;i++){ times[i]=tmin+i*dt; addtr[i]=0.; }

	if (!dt) CSMUSTGETPARFLOAT("dt", &dt);

	if(!parObj.getparstring("type",&type)){ type="nmo"; typecode=0; }
	if(!strcmp(type,"lmo"))typecode=1;

	if(!parObj.getparfloat("t0",&t0))t0=1.;
	if(!parObj.getparfloat("vel",&vel))vel=3000.;
	if(!parObj.getparfloat("amp",&amp))amp=1.;

	while(nsegy){
		if(typecode==1) /* lmo */ tx=(t0+abs(tr.offset)/vel);   
		else /* nmo */ tx=sqrt( SQ(t0)+SQ((float)(tr.offset/vel)) ); 

		ints8r( 1, dt, tx, &amp, 0.0, 0.0, nt, times, addtr );

		for(i=0;i<nt;i++){
			tr.data[i]+=addtr[i];
			addtr[i]=0;
		}

		su2cs->putTrace(&tr);
		nsegy=cs2su->getTrace(&tr);
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
