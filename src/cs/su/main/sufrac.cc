/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUFRAC: $Revision: 1.29 $ ; $Date: 2011/12/21 23:25:42 $	*/

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
std::string sdoc_sufrac =
"									"
" SUFRAC -- take general (fractional) time derivative or integral of	"
"	    data, plus a phase shift.  Input is CAUSAL time-indexed	"
"	    or depth-indexed data.					"
"									"
" sufrac power= [optional parameters] <indata >outdata 			"
"									"
" Optional parameters:							"
"	power=0		exponent of (-i*omega)	 			"
"			=0  ==> phase shift only			"
"			>0  ==> differentiation				"
"			<0  ==> integration				"
"									"
"	sign=-1			sign in front of i * omega		"
"	dt=(from header)	time sample interval (in seconds)	"
"	phasefac=0		phase shift by phase=phasefac*PI	"
"	verbose=0		=1 for advisory messages		"
"									"
" Examples:								"
"  preprocess to correct 3D data for 2.5D migration			"
"         sufrac < sudata power=.5 sign=1 | ...				"
"  preprocess to correct susynlv, susynvxz, etc. (2D data) for 2D migration"
"         sufrac < sudata phasefac=.25 | ...				"
" The filter is applied in frequency domain.				"
" if dt is not set in header, then dt is mandatory			"
"									"
" Algorithm:								"
"		g(t) = Re[INVFTT{ ( (sign) iw)^power FFT(f)}]		"
" Caveat:								"
" Large amplitude errors will result if the data set has too few points."
"									"
" Good numerical integration routine needed!				"
" For example, see Gnu Scientific Library.				"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sufrac {


/* Credits:
 *	CWP: Chris Liner, Jack K. Cohen, Dave Hale (pfas)
 *      CWP: Zhenyue Liu and John Stockwell added phase shift option
 *	CENPET: Werner M. Heigl - added well log support
 *
 * Trace header fields accessed: ns, dt, trid, d1
*/
/**************** end self doc ***********************************/


#define	I		cmplx(0.0, 1.0)
#define	PIBY2		0.5 * PI
#define TWOPI		2.0 * PI
#define LOOKFAC		2	/* Look ahead factor for npfao	  */
#define PFA_MAX		720720	/* Largest allowed nfft	          */

segy tr;

void* main_sufrac( void* args )
{
	float power;		/* power of i omega applied to data	*/
	float amp;		/* amplitude associated with the power	*/
	float arg;		/* argument of power 			*/
	float phasefac;		/* phase factor	 			*/
	float phase;		/* phase shift = phasefac*PI		*/
	complex exparg;		/* cwp_cexp(I arg)				*/
	register float *rt;	/* real trace				*/
	register complex *ct;	/* complex transformed trace		*/
	complex *filt;		/* complex power	 		*/
	int nt;			/* number of points on input trace	*/
	size_t ntsize;		/* nt in bytes				*/
	float dt;		/* sample spacing (secs) on input trace	*/
	float omega;		/* circular frequency			*/
	float domega;		/* circular frequency spacing (from dt)	*/
	float sign;		/* sign in front of i*omega default -1	*/
	int nfft;		/* number of points in nfft		*/
        int nf;                 /* number of frequencies (incl Nyq)     */
	float onfft;		/* 1 / nfft				*/
	int verbose;		/* flag to get advisory messages	*/
	size_t nzeros;		/* number of padded zeroes in bytes	*/
	cwp_Bool seismic;	/* is this seismic data?		*/
	
	
	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sufrac );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Set parameters */
	if (!parObj.getparint("verbose", &verbose))	  verbose  =  0;
	if (!parObj.getparfloat("power", &power))	  power    =  0.0;
	if (!parObj.getparfloat("sign", &sign))	  sign     = -1.0;
	if (!parObj.getparfloat("phasefac", &phasefac))  phasefac =  0.0;
	phase = phasefac * PI;


	/* Get info from first trace */
	if (!cs2su->getTrace(&tr))	throw cseis_geolib::csException("can't get first trace");
	seismic = ISSEISMIC(tr.trid);
	if (seismic) {
		if (verbose)	warn("input is seismic data, trid=%d",tr.trid);
		dt = ((double) tr.dt)/1000000.0;
	}
	else {
		if (verbose)	warn("input is not seismic data, trid=%d",tr.trid);
		dt = tr.d1;
        }
        if (!dt)	throw cseis_geolib::csException("dt or d1 field is zero and not parObj.getparred");
	nt = tr.ns;
	ntsize = nt * FSIZE;


	/* Set up for fft */
	nfft = npfaro(nt, LOOKFAC * nt);
	if (nfft >= SU_NFLTS || nfft >= PFA_MAX)
		throw cseis_geolib::csException("Padded nt=%d -- too big", nfft);

        nf = nfft/2 + 1;
        onfft = 1.0 / nfft;
	nzeros = (nfft - nt) * FSIZE;
	domega = TWOPI * onfft / dt;


	/* Allocate fft arrays */
	rt   = ealloc1float(nfft);
	ct   = ealloc1complex(nf);
	filt = ealloc1complex(nf);


	/* Set up args for complex power evaluation */
	arg = sign * PIBY2 * power + phase;
	exparg = cwp_cexp(crmul(I, arg));


	/* Evaluate complex power, put inverse fft scale in */
	{
		register int i;
		for (i = 0 ; i < nf; ++i) {

			omega = i * domega;

			/* kludge to handle omega=0 case for power < 0 */
			if (power < 0 && i == 0) omega = FLT_MAX;

			/* calculate filter */
			amp = pow(omega, power) * onfft;

			filt[i] = crmul(exparg, amp);
		}
	}
		

	/* Loop over traces */
	do {
		/* Load trace into rt (zero-padded) */
		memcpy( (void *) rt, (const void *) tr.data, ntsize);
		memset((void *) (rt + nt),0, nzeros);

		/* FFT */
		pfarc(1, nfft, rt, ct);


		/* Apply filter */
		{ register int i;
		for (i = 0; i < nf; ++i)  ct[i] = cmul(ct[i], filt[i]);
		}

		/* Invert */
		pfacr(-1, nfft, ct, rt);


		/* Load traces back in, recall filter had nfft factor */
		{ register int i;
		for (i = 0; i < nt; ++i)  tr.data[i] = rt[i];
		}


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
