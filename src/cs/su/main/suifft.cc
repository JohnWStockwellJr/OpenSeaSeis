/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUIFFT: $Revision: 1.15 $ ; $Date: 2011/11/16 23:35:04 $	*/

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
std::string sdoc_suifft =
" 								"
" SUIFFT - fft complex frequency traces to real time traces	"
" 								"
" suiftt <stdin >sdout sign=-1					"
" 								"
" Required parameters:						"
" 	none							"
" 								"
" Optional parameter:						"
" 	sign=-1		sign in exponent of inverse fft		"
" 								"
" Output traces are normalized by 1/N where N is the fft size.	"
" 								"
" Note: sufft | suifft is not quite a no-op since the trace	"
" 	length will usually be longer due to fft padding.	"
" 								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suifft {


/* Credits:
 *
 *	CWP: Shuki, Chris, Jack
 *
 * Trace header fields accessed: ns, trid
 * Trace header fields modified: ns, trid
 */
/**************** end self doc ***********************************/


#define PFA_MAX	720720		/* Largest allowed fft	*/

segy tr;

void* main_suifft( void* args )
{
	register complex *ct;	/* complex input trace			*/
	register float *rt;	/* real output trace			*/
	int nfft;		/* fft size 				*/
	int nf;			/* number of frequencies		*/
	int sign;		/* sign in exponent of transform	*/
	float onfft;		/* 1.0/nfft				*/
	float newd1;		/* reset time domain sampling		*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suifft );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get info from first trace */ 
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	if (tr.trid != FUNPACKNYQ)
		throw cseis_geolib::csException("input not complex freq data, trid=%d", tr.trid);
	nfft = tr.ns - 2; /* see sufft */
	nf = nfft/2 + 1;
	onfft = 1.0/nfft;
	if(tr.d1) {
	    newd1 = onfft/tr.d1;
	} else {
	   if(tr.dt) newd1 = (float) (((double) tr.dt)/1000000.0);
	   else newd1 = 0.0f;
	}



	/* Set sign in exponent of transform */
	if (!parObj.getparint   ("sign", &sign)) sign = -1;
	if (sign != 1 && sign != -1)  throw cseis_geolib::csException("sign = %d must be 1 or -1", sign);


        parObj.checkpars();


	/* Allocate fft arrays */
	ct   = ealloc1complex(nf);
	rt   = ealloc1float(nfft);


	/* Main loop over traces */
	do {
		register int i;

		/* Load traces into ct (pfa fills in negative freqs) */
		for (i = 0; i < nf; ++i) {
			ct[i].r = tr.data[2*i];
			ct[i].i = tr.data[2*i+1];
		}


		/* Inverse FFT */
		pfacr(sign, nfft, ct, rt);


		/* Load back and scale for inverse fft */
		for (i = 0; i < nfft; i++)  tr.data[i] = rt[i] * onfft;

		tr.trid = TREAL;
		tr.ns = nfft;
		tr.f1 = 0.0f;
		tr.d1 = newd1;

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
