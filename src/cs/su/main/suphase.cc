/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUPHASE: $Revision: 1.4 $ ; $Date: 2011/11/16 18:03:07 $		*/

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
std::string sdoc_suphase =
" 									"
" SUPHASE - PHASE manipulation by linear transformation			"
" 									"
"  suphase  <stdin >sdout      						"
" 									"
" Required parameters:							"
" none									"
" Optional parameters:							"
" a=90			constant phase shift              		"
" b=180/PI              linear phase shift				"
" c=0.0			phase = a +b*(old_phase)+c*f;			"
" 									"
" Notes: 								"
" A program that allows the user to experiment with changes in the phase"
" spectrum of a signal.							"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suphase {


/**************** end self doc ********************************/

/* definitions used internally */
#define AMPSP(c) sqrt(c.r*c.r+c.i*c.i)
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

/* segy trace */
segy tr;

void* main_suphase( void* args )
{
	float *rt=NULL;		/* real trace			*/
	float *amp=NULL;	/* amplitude spectra		*/
	float *ph=NULL;		/* phase			*/
	register complex *ct=NULL;	/* complex time trace	*/

	int nt;			/* number of points on input trace	*/
	int nfft;		/* transform length			*/
	int nf;			/* number of frequencies in transform	*/

	float dt;		/* sampling interval in secs		*/
	float d1;		/* output sample interval in Hz		*/
	int count=0;		/* counter				*/

	/* linear phase function */
	float a;		/* bias (intercept) of new phase	*/
	float b;		/* slope of linear phase function	*/
	float c;		/* new phase value			*/
	float onfft;		/* 1/nfft				*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suphase );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get info from first trace */ 
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;

	/* get parameters */
	/* dt is used only to set output header value d1 */
	if (!parObj.getparfloat("dt", &dt))	dt = ((double) tr.dt)/1000000.0;
	if (!dt) {
		dt = .004;
		warn("dt not set, assumed to be .002");
	}

	/* linear phase paramter values */
	if (!parObj.getparfloat("a", &a)) a = 0;
	if (!parObj.getparfloat("b", &b)) b = 180/PI;
	if (!parObj.getparfloat("c", &c)) c = 0.0;

	a *= PI/180.0;
	b *= PI/180.0;
	
	/* Set up pfa fft */
	nfft = npfaro(nt, LOOKFAC * nt);
	if (nfft >= SU_NFLTS || nfft >= PFA_MAX)  
		throw cseis_geolib::csException("Padded nt=%d--too big", nfft);
	d1 = 1.0/(nfft*dt);
	nf = nfft/2 + 1;
        onfft = 1.0/nfft;
	 
	parObj.checkpars();

	/* Allocate space */
	rt = ealloc1float(nfft);
	ct = ealloc1complex(nf);
	amp = ealloc1float(nf);
	ph = ealloc1float(nf);

	/* Main loop over traces */
	count=0;
	do {
		register int i;
		
		/* Load trace into rt (zero-padded) */
		memcpy((void *) rt, (const void *) &tr.data, nt*FSIZE);
		memset((void *) (rt + nt), (int) '\0', (nfft-nt)*FSIZE);
		
		/* FFT */
		pfarc(1, nfft, rt, ct);
		for (i = 0; i < nf; ++i) {
			amp[i] = AMPSP(ct[i]);
			ph[i]  = a+b*atan2(ct[i].i,ct[i].r)+c*i;
		}
		for (i = 0; i < nf; ++i) {
			ct[i].r = amp[i]*cos(ph[i]);
			ct[i].i = amp[i]*sin(ph[i]);
		}
		pfacr(-1,nfft,ct,rt);
		for (i = 0; i < nt; ++i) rt[i]*=onfft;
		memcpy((void *) tr.data, (const void *) rt, nt*FSIZE);
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
