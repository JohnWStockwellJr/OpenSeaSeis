/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.		       */

/* SUAMP: $Revision: 1.22 $ ; $Date: 2011/11/23 22:27:07 $		*/

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
std::string sdoc_suamp =
" 									"
" SUAMP - output amp, phase, real or imag trace from			"
" 	(frequency, x) domain data					"
" 									"
" suamp <stdin >stdout mode=amp						"
" 									"
" Required parameters:							"
" none									"
" Optional parameter:							"
" mode=amp	output flag		 				"
" 		=amp	output amplitude traces				"
" 	       	=phase	output phase traces				"
" 	       	=ouphase output unwrapped phase traces (oppenheim)	"
" 	       	=suphase output unwrapped phase traces (simple)		"
" 	       	=real	output real parts				"
" 	     	=imag	output imag parts	 			"
" jack=0	=1  divide value at zero frequency by 2   		"
"		(operative only for mode=amp)				"
" .... phase unwrapping options	..... 					"
" unwrap=1	 |dphase| > pi/unwrap constitutes a phase wrapping	"
"			(operative only for mode=suphase)		"
" trend=1 	 remove linear trend from the unwrapped phase		"
" zeromean=0 	 assume phase(0)=0.0, else assume phase is zero mean	"
" Notes:								"
" 	The trace returned is half length from 0 to Nyquist. 		"
" 									"
" Example:								"
" 	sufft <data | suamp >amp_traces					"
" Example: 								"
"	sufft < data > complex_traces					"
" 	 suamp < complex_traces mode=real > real_traces			"
" 	 suamp < complex_traces mode=imag > imag_traces			"
"  									"
" Note: the inverse of the above operation is: 				"
"       suop2 real_traces imag_traces op=zipper > complex_traces	"
"  									"
" Note: Explanation of jack=1 						"
" The amplitude spectrum is the modulus of the complex output of	"
" the fft. f(0) is thus the average over the range of integration	"
" of the transform. For causal functions, or equivalently, half		"
" transforms, f(0) is 1/2 of the average over the full range.		"
" Most oscillatory functions encountered in wave applications are	"
" zero mean, so this is usually not an issue.				"
"  									"
" Note: Phase unwrapping: 						"
"  									"
" The mode=ouphase uses the phase unwrapping method of Oppenheim and	"
" Schaffer, 1975. 							"
" The mode=suphase generates unwrapped phase assuming that jumps	"
" in phase larger than pi/unwrap constitute a phase wrapping.		"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suamp {


/* Credits:
 *	CWP: Shuki Ronen, Jack K. Cohen c.1986
 *
 * Notes:
 *	If efficiency becomes important consider inverting main loop
 *      and repeating extraction code within the branches of the switch.
 *
 * Trace header fields accessed: ns, trid
 * Trace header fields modified: ns, trid
 */
/**************** end self doc ***********************************/


#define	REAL	1
#define	IMAG	2
#define	AMP	3
#define	ARG	4
#define	SUPHASE	5
#define	OUPHASE	6

segy tr;

void* main_suamp( void* args )
{
	cwp_String mode;	/* display: real, imag, amp, arg	*/
	int imode=AMP;		/* integer abbrev. for mode in switch	*/
	int nfby2;		/* nf/2					*/
	register float *xr;	/* real part of trace			*/
	register float *xi;	/* imaginary part of trace		*/
	int jack=0;		/* flag for special treatment of zero omega */
	float unwrap;		/* PI/unwrap = min dphase assumed to wrap */

	int zeromean;		/* =0 assume phase(0)=0.0 ; =1  zero mean*/
	int verbose=0;		/* =1 chatty ; =0 silent		*/
	int trend ;		/* remove linear trend in phase unwrapping */

	float df;		/* frequency sampling interval		*/
	int nf;			/* number of samples on input trace	*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suamp );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get information from first trace */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	if (tr.trid != FUNPACKNYQ)
		throw cseis_geolib::csException("input not complex freq data, trid=%d", tr.trid);

	nf = tr.ns; /* always even--includes 1 or 2 "wasted" slots */
	nfby2 = nf/2;
	if (!parObj.getparint("verbose",&verbose))	verbose=0 ;

	df=tr.d1;
	if (!df) {
		float dt = ((double) tr.dt)/1000000.0;
		if (!dt) {
			dt = .004;
			if (verbose) warn("dt or d1 not set, assumed to be .004");
		}
		df = 1.0/nf*dt;
	}


	/* Get mode; note that imode is initialized to AMP */
	if (!parObj.getparstring("mode", &mode))	mode = "amp";
	if (!parObj.getparint("jack",&jack))	jack = 0;


	if      (STREQ(mode, "phase")) imode = ARG;
	else if (STREQ(mode, "ouphase"))  imode = OUPHASE;
	else if (STREQ(mode, "suphase"))  imode = SUPHASE;
	else if (STREQ(mode, "real"))  imode = REAL;
	else if (STREQ(mode, "imag"))  imode = IMAG;
	else if (!STREQ(mode, "amp"))
		throw cseis_geolib::csException("unknown mode=\"%s\", see self-doc", mode);


	if(imode==OUPHASE || imode==SUPHASE) {
		if (!parObj.getparint("trend",&trend))	trend = 1;
		if (!parObj.getparint("zeromean",&zeromean))	zeromean = 0;
		if (!parObj.getparfloat("unwrap",&unwrap))	unwrap = 1;
	}
		

	parObj.checkpars();

	/* Allocate arrays for real and imaginary parts */
	xr = ealloc1float(nfby2);
	xi = ealloc1float(nfby2);

	/* Main loop over traces */
	do {
		register int i;

		/* Separate complex trace into real and imag parts */
		for (i = 0; i < nfby2; ++i) {
			xr[i] = tr.data[2*i];
			xi[i] = tr.data[2*i+1];
		}

		/* Compute the desired half-length trace */
		switch(imode) {
		case REAL:
			for (i = 0; i < nfby2; ++i) {
				tr.data[i] = xr[i];
			}
			tr.trid = REALPART;
		break;
		case IMAG:
			for (i = 0; i < nfby2; ++i) {
				tr.data[i] = xi[i];
			}
			tr.trid = IMAGPART;
		break;
		case AMP:
		{
	 		register float re, im;
	
			re = xr[0];
			im = xi[0];
			if (jack) {
				tr.data[0] = (float) sqrt (re * re + im * im) / 2.0;
			} else {
				tr.data[0] = (float) sqrt (re * re + im * im);
			}
				for (i = 1; i < nfby2; ++i) {
					re = xr[i];
					im = xi[i];
					tr.data[i] = (float) sqrt (re * re + im * im);
			}
			tr.trid = AMPLITUDE;
		}
		break;
		case ARG:
			for (i = 0; i < nfby2; ++i) {
				float re = xr[i];
				float im = xi[i];

				if (re*re+im*im) {
					tr.data[i] = atan2(im, re);
				} else {
					tr.data[i] = 0.0;
				}
			}
			tr.trid = PHASE;
		break;
		case SUPHASE:
		{
			float *phase=NULL;
			
			/* allocate space for the phase */
			phase = alloc1float(nfby2);

			for (i = 0; i < nfby2; ++i) {
				float re = xr[i];
				float im = xi[i];

				if (re*re+im*im) {
					phase[i] = atan2(im, re);
				} else {
					phase[i] = 0.0;
				}
			}

			/* unwrap the phase */
			if (unwrap)
				simple_unwrap_phase(nfby2, trend, zeromean, unwrap, phase);

			/* write unwrapped phase */
			for ( i = 0; i < nfby2; ++i) tr.data[i] = phase[i];

			tr.trid = PHASE;
		}
		break;
		case OUPHASE:
		{
			float *phase=NULL;
			
			/* allocate space */
			phase = ealloc1float(nfby2);
			
			memset((void *) phase,  0, nfby2*FSIZE);
			/* unwrap the phase */
			if (unwrap)
				oppenheim_unwrap_phase(nfby2, trend, zeromean, df, xr, xi, phase);

			/* write unwrapped phase */
			for ( i = 0; i < nfby2; ++i) tr.data[i] = phase[i];

			tr.trid = PHASE;
		}
		break;
		default:
			throw cseis_geolib::csException("mysterious mode=\"%s\"", mode);
		}

		/* Output the half-length trace */
		tr.ns = nfby2;
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
