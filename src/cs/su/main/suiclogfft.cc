/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUICLOGFFT: $Revision: 1.8 $ ; $Date: 2011/11/16 23:35:04 $	*/

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
std::string sdoc_suiclogfft =
" 								"
" SUICLOGFFT - fft of complex log frequency traces to real time traces"
" 								"
" suiclogftt <stdin >sdout sign=-1				"
" 								"
" Required parameters:						"
" 	none							"
" 								"
" Optional parameter:						"
" 	sign=-1		sign in exponent of inverse fft		"
"	sym=0		=1 center  output 			"
" Output traces are normalized by 1/N where N is the fft size.	"
" 								"
" Note:								"
" Nominally this is the inverse to the complex log fft, but	"
" suclogfft | suiclogfft is not quite a no-op since the trace	"
" 	length will usually be longer due to fft padding.	"
" 								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suiclogfft {


/* Credits:
 * 
 *   CWP: John Stockwell, Dec 2010 based on
 *     suifft.c by:
 *	CWP: Shuki Ronen, Chris Liner, Jack K. Cohen,  c. 1989
 *
 * Trace header fields accessed: ns, trid
 * Trace header fields modified: ns, trid
 */
/**************** end self doc ***********************************/


#define PFA_MAX	720720		/* Largest allowed fft	*/

segy tr;

void* main_suiclogfft( void* args )
{
	register complex *ct;	/* complex input trace			*/
	register float *rt;	/* real output trace			*/
	int nfft;		/* fft size 				*/
	int nf;			/* number of frequencies		*/
	int sign;		/* sign in exponent of transform	*/
	float onfft;		/* 1.0/nfft				*/
	float newd1;		/* reset time domain sampling		*/

	int sym=0;		/* =1 symmetric output			*/

	int nfftby2;		/* nfft/2				*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suiclogfft );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get info from first trace */ 
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	if (tr.trid != FUNPACKNYQ)
		throw cseis_geolib::csException("input not complex freq data, trid=%d", tr.trid);
	nfft = tr.ns - 2; /* see sufft */
	nfftby2 = nfft/2;
	nf = nfft/2 + 1;
	onfft = 1.0/nfft;
	if(tr.d1) {
	    newd1 = onfft/tr.d1;
	} else {
	   if(tr.dt) newd1 = (float) (((double) tr.dt)/1000000.0);
	   else newd1 = 0.0f;
	}

	/* Set sign in exponent of transform */
	if (!parObj.getparint   ("sym", &sym)) sym = 0;
	if (!parObj.getparint   ("sign", &sign)) sign = -1;

        parObj.checkpars();

	if (sign != 1 && sign != -1)  throw cseis_geolib::csException("sign = %d must be 1 or -1", sign);


	/* Allocate fft arrays */
	ct   = ealloc1complex(nf);
	rt   = ealloc1float(nfft);


	/* Main loop over traces */
	do {
		register int i;

		/* zero out array */
		memset( (void *) ct, 0, nf*sizeof(complex));
		memset( (void *) rt, 0, nfft*FSIZE);
		
		/* Load traces into ct (pfa fills in negative freqs) */
		/* exponentiate prior to inverse transform */
		for (i = 0; i < nf; ++i) {
			if (tr.data[2*i]) {
				ct[i].r = (float) exp(tr.data[2*i]) * cos(tr.data[2*i+1]);
				ct[i].i = (float) exp(tr.data[2*i]) * sin(tr.data[2*i+1]);
			} else {

				ct[i].r = 0.0;
				ct[i].i = 0.0;
			}
		}


		/* Inverse FFT */
		pfacr(sign, nfft, ct, rt);

		/* Load back and scale for inverse fft */
		if (!sym) {
			for (i = 0; i < nfft; ++i)  
				tr.data[i] = rt[i] * onfft;
		} else {/* symmetric output */
			for(i = nfftby2  ; i < nfft ; ++i )
				tr.data[i] = rt[i-nfftby2] * onfft;

			for(i = 0; i < nf; ++i) 
				tr.data[i] = rt[i+nf] * onfft;
		}

		tr.trid = TREAL;
		tr.ns = nfft;
		tr.f1 = 0.0f;
		tr.d1 = newd1;
		if (sym) tr.delrt=-(tr.delrt + (nfft*newd1*1000.0)/2.0);

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
