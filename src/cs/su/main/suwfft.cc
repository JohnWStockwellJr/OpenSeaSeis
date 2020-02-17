/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUWFFT: $Revision: 1.4 $ ; $Date: 2011/11/16 23:35:04 $		*/

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
std::string sdoc_suwfft =
" 									"
" SUWFFT - Weighted amplitude FFT with spectrum flattening 0->Nyquist	"
" 									"
" suwfft <stdin | suifft >sdout 					"
" 									"
" Required parameters:							"
" none									"
" 									"
" Optional parameters:							"
" w0=0.75		weight for AmpSpectrum[f-df]			"
" w1=1.00		weight for AmpSpectrum[f].. center value	"
" w2=0.75		weight for AmpSpectrum[f+df]			"
" 									"
" Notes: 								"
" 1. output format is same as sufft					"
" 2. suwfft | suifft is not quite a no-op since the trace		"
"    length will usually be longer due to fft padding.			"
" 3. using w0=0 w1=1 w2=0  gives truly flat spectrum, for other	        "
"    weight choices the spectrum retains some of its original topograpy "
" 									"
" Examples: 								"
" 1. boost data bandwidth to 10-90 Hz					"
"     suwfft < data.su | suifft | sufilter f=5,8,90,100 | suximage 	"
" 1. view amplitude spectrum after flattening				"
"     suwfft < data.su | suamp | suximage 				"
" 									"
" Caveat: The process of cascading the forward and inverse Fourier	"
"  transforms may result in the output trace length being greater than 	"
"  the input trace length owing to zero padding. The user may wish to	"
"  apply suwind to return the number of samples per trace to the original"
"  value:  Here NS is the number of samples per trace on the original data"
"  			... | suwind itmax=NS | ... 			"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suwfft {


/* Credits:
 *
 *	UHouston: Chris Liner 
 *
 * Note: Concept from UTulsa PhD thesis of Bassel Al-Moughraby
 *
 * Trace header fields accessed: ns, dt
 * Trace header fields modified: ns, d1, f1, trid
 */
/**************** end self doc ***********************************/



#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

segy tr;

void* main_suwfft( void* args )
{
	register float *rt;	/* real trace				*/
	register complex *ct;	/* complex transformed trace		*/
	int nt;			/* number of points on input trace	*/
	int nfft;		/* transform length			*/
	int nf;			/* number of frequencies		*/
	int sign;		/* sign in exponent of transform	*/
	int verbose;		/* flag to get advisory messages	*/
	float dt;		/* sampling interval in secs		*/
	float d1;		/* output sample interval in Hz		*/
	cwp_Bool seismic;	/* is this seismic data? */
	float c;		/* multiplier				*/
	float w0, w1, w2;	/* weights				*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suwfft );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	if (!parObj.getparint("verbose", &verbose))	verbose=1;

	/* Get info from first trace */ 
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;

	/* check for seismic or well log data */
	seismic = ISSEISMIC(tr.trid);		
	if (seismic) {
		if (verbose)	warn("input is seismic data, trid=%d",tr.trid);
		dt = ((double) tr.dt)/1000000.0;
	}
	else {
		if (verbose)	warn("input is not seismic data, trid=%d",tr.trid);
		dt = tr.d1;
        }
	if (!dt) {
		dt = .004;
		if (verbose) warn("dt or d1 not set, assumed to be .004");
	}

	/* Set up pfa fft */
	nfft = npfaro(nt, LOOKFAC * nt);
	if (nfft >= SU_NFLTS || nfft >= PFA_MAX)  throw cseis_geolib::csException("Padded nt=%d--too big", nfft);
	nf = nfft/2 + 1;
	d1 = 1.0/(nfft*dt);

	if (!parObj.getparint("sign", &sign)) sign = 1;
	if (sign != 1 && sign != -1)   throw cseis_geolib::csException("sign = %d must be 1 or -1", sign);

	/* get weights */
	if (!parObj.getparfloat("w0",&w0)) w0 = 0.75;
	if (!parObj.getparfloat("w1",&w1)) w1 = 1.00;
	if (!parObj.getparfloat("w2",&w2)) w2 = 0.75;


        parObj.checkpars();

	rt = ealloc1float(nfft);
	ct = ealloc1complex(nf);

	/* If dt not set, issue advisory on frequency step d1 */
	if (dt && verbose)  warn("d1=%f", 1.0/(nfft*dt));

	/* Main loop over traces */
	do {
		register int i;

		/* Load trace into rt (zero-padded) */
		memcpy((void *) rt, (const void *) tr.data, nt*FSIZE);
		memset((void *) (rt + nt), (int) '\0', (nfft-nt)*FSIZE);

		/* FFT */
		pfarc(sign, nfft, rt, ct);

		/* Store values */
		for (i = 0; i < nf; ++i) {
			c =w0*rcabs(ct[i-1])+w1*rcabs(ct[i])+w2*rcabs(ct[i+1]);
			if (i==0 || i==nf) {
				tr.data[2*i]   = ct[i].r / rcabs(ct[i]);
				tr.data[2*i+1] = ct[i].i / rcabs(ct[i]);
			} else {
				tr.data[2*i]   = ct[i].r / c;
				tr.data[2*i+1] = ct[i].i / c;
			}
		}

		/* Set header values--npfaro makes nfft even */
		tr.ns = 2 * nf;
		tr.trid = FUNPACKNYQ;
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

} // END namespace
