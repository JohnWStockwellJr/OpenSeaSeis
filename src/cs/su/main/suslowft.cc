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

/*********************** self documentation **********************/
std::string sdoc_suslowft =
" 								"
" SUSLOWFT - Fourier Transforms by a (SLOW) DFT algorithm (Not an FFT)"
" 								"
" suslowft <stdin >sdout sign=1 				"
" 								"
" Required parameters:						"
" 	none							"
" 								"
" Optional parameters:						"
" 	sign=1			sign in exponent of fft		"
" 	dt=from header		sampling interval		"
" 								"
" Trace header fields accessed: ns, dt				"
" Trace header fields modified: ns, dt, trid			"
" 								"
" Notes: To facilitate further processing, the sampling interval"
"       in frequency and first frequency (0) are set in the	"
"	output header.						"
" Warning: This program is *not* fft based. Use only for demo 	"
" 	   purposes, *not* for large data processing.		"
" 								"
" 	No check is made that the data are real time traces!	"
" suslowft | suslowift is not quite a no-op since the trace     "
" length will usually be longer due to fft padding.             "
"                                                                       "
" Caveats:                                                              "
" No check is made that the data IS real time traces!                   "
"                                                                       "
" Output is type complex. To view amplitude, phase or real, imaginary   "
" parts, use    suamp                                                   "
"                                                                       "
" Examples:                                                             "
" suslowft < stdin | suamp mode=amp | ....                                 "
" suslowft < stdin | suamp mode=phase | ....                               "
" suslowft < stdin | suamp mode=real | ....                                "
" suslowft < stdin | suamp mode=imag | ....                                "
"                                                                       "
" 								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suslowft {


/* Credits:
 *
 *	CWP: Shuki, Chris, Jack
 *
 * Note: leave dt set for later inversion
 *
 */
/**************** end self doc ***********************************/

/* Prototype of function used internally */
void dftcc (int sign, int nsamp, complex *cz);
void dftrc (int sign, int nsamp, float *re, complex *out);

segy tr;

void* main_suslowft( void* args )
{
	register float *rt;	/* real trace				*/
	register complex *ct;	/* complex transformed trace		*/
	int nt;			/* number of points on input trace	*/
	int sign;		/* sign in exponent of transform	*/
	float nfft;		/* fft size				*/
	float dt;		/* sampling interval in secs		*/
	float d1;		/* output sample interval in Hz		*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suslowft );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get info from first trace */ 
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;
	nfft=nt;

	/* dt is used only to set output header value d1 */
	if (!parObj.getparfloat("dt", &dt)) dt = (float) tr.dt/1000000.0;
	if (!dt) {
		dt = .004;
		warn("dt not set, assumed to be .004");
	}

	d1 = 1.0/(nt*dt);

	if (!parObj.getparint("sign", &sign)) sign = 1;
	if (sign != 1 && sign != -1)   throw cseis_geolib::csException("sign = %d must be 1 or -1", sign);


        parObj.checkpars();

	/* allocate space */
	rt = ealloc1float(nfft);
	ct = ealloc1complex(nfft);

	/* Echo frequency step size (Hz) to par file */
	if (dt)  warn("d1=%f f1=0.0", 1.0/(nt*dt));

	/* Main loop over traces */
	do {
		register int i;

		/* Load trace into rt  */
                memcpy((void *) rt, (const void *) tr.data, nt*FSIZE);


		/* FFT */
		dftrc(sign, nfft, rt, ct);

		/* Store values */
		for (i = 0; i < nfft; ++i) {
			tr.data[2*i]   = ct[i].r;
			tr.data[2*i+1] = ct[i].i;
		}

		/* Set header values */
		tr.trid = FUNPACKNYQ;
		tr.ns = nt;
		tr.d1 = d1;
		tr.f1 = 0.0;

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

void dftcc (int sign, int nsamp, complex *cz)
/*
 * Digital Fourier Transform --- (Not an FFT) don't use if speed is an issue.
 *
 * Author: John Stockwell, 26 Feb 1992, based on the well known algorithm
 *
 * 	See: discussion in Numerical Recipes in C, p.406 or other equivalent
*/
{
	int ki,ni;				/* looping variables 	*/
	float pibynsamp = (2.* PI/nsamp);	/* exponent multiplier 	*/
	float tmpcos,tmpsin;	/* temporary variables 	*/
	float *tmpim=NULL,*tmpre=NULL;		/* temporary array pointers*/ 
	
	/* allocate space for temp arrays */
	tmpim = ealloc1float(nsamp);
	tmpre = ealloc1float(nsamp);
	
	for (ni=0; ni<nsamp; ++ni) {
		tmpim[ni] = 0.0;
		tmpre[ni] = 0.0;
		cz[ni] = cmplx(0.0,0.0);
	}

	/* do the dftcc */
	for(ni=0; ni<nsamp; ni++) {
		register float retemp = 0.0;
		register float imtemp = 0.0;

		for(ki=0; ki<nsamp; ki++) {

			tmpcos = cos(sign*pibynsamp * ki * ni);
			tmpsin = sin(sign*pibynsamp * ki * ni);

			retemp = retemp + (cz[ki].r + cz[ki].i)*(tmpcos - tmpsin);
			imtemp = imtemp + (cz[ki].r*tmpsin + cz[ki].i*tmpcos);
		}
		tmpre[ni] = retemp;
		tmpim[ni] = imtemp;
	}

	/* prepare new values of re[] and im[] */
	for(ni=0; ni<nsamp; ni++) {
		cz[ni].r=tmpre[ni];
		cz[ni].i=tmpim[ni];
	}
}
#include "cwp.h"
void dftrc (int sign, int nsamp, float *re, complex *out)
/*
 * Discrete Fourier Transform --- (Not an FFT) don't use if speed is an issue.
 *
 * Author: John Stockwell, 26 Feb 1992, based on the well known algorithm
 *
 * 	See: Numerical Recipes in C, p.406 or other 
 */

{
	int ki,ni;		/* looping variables 	*/
	float freq,phase;	/* argument to cos and sin */
	float twopibynsamp=0.0;	/* 2*PI/nsamp */
	float *tmpre;
	float *tmpim;

	/* allocate space for temporary arrays */
	tmpre = ealloc1float(nsamp);
	tmpim = ealloc1float(nsamp);

	/* zero out arrays */
	for (ni=0; ni<nsamp; ++ni) {
		tmpim[ni] = 0.0;
		tmpre[ni] = 0.0;
		out[ni] = cmplx(0.0,0.0);
	}
	
	twopibynsamp = sign*2*PI/nsamp;
	freq = phase = 0.0;

	/* do the dft assuming that the input are real data only */
	for(ni=0; ni<nsamp; ++ni) {
		register float retemp = 0.0;
		register float imtemp = 0.0;
		
		
		for(ki=0; ki<nsamp; ++ki) {
			register float tmpcos=0.0;
			register float tmpsin=0.0;
			
			phase = freq*ki;

			tmpcos = cos(phase);
			tmpsin = sin(phase);

/* 			retemp = retemp + re[ki]*(tmpcos - tmpsin);
			imtemp = imtemp + (re[ki]*tmpsin);
 */
			retemp += re[ki]*tmpcos;
			imtemp += re[ki]*tmpsin;
			
			/*phase += freq;*/
		}
		tmpre[ni] = retemp;
		tmpim[ni] = imtemp;
		
		phase = 0.0;
		freq += twopibynsamp; /* saves a float mult */
	}

	/* prepare new values of re[] and im[] */
	for(ni=0; ni<nsamp; ni++) {
		out[ni].r = tmpre[ni];
		out[ni].i = tmpim[ni];
	}

}

} // END namespace
