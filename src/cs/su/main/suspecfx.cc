/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSPECFX: $Revision: 1.19 $ ; $Date: 2011/11/16 23:35:04 $		*/

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
std::string sdoc_suspecfx =
" 								"
" SUSPECFX - Fourier SPECtrum (T -> F) of traces 		"
" 								"
" suspecfx <infile >outfile 					"
" 								"
" Note: To facilitate further processing, the sampling interval	"
"       in frequency and first frequency (0) are set in the	"
"	output header.						"
" 								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suspecfx {


/* Credits:
 *
 *	CWP: Dave (algorithm), Jack (reformatting for SU)
 *
 * Trace header fields accessed: ns, dt
 * Trace header fields modified: ns, dt, trid, d1, f1
 */
/**************** end self doc ***********************************/


#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

segy tr;

void* main_suspecfx( void* args )
{
	register float *rt;	/* real trace				*/
	register complex *ct;	/* complex transformed trace		*/
	int nt;			/* number of points on input trace	*/
	int nfft;		/* number of points on output trace	*/
	int nfby2p1;		/* nfft/2 + 1				*/
	float dt;		/* sample interval in secs		*/
	float d1;		/* output sample interval in Hz		*/
	int ntr=0;		/* number of traces			*/
	register int i;		/* counter				*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suspecfx );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get info from first trace */ 
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;

	/* dt is used only to set output header value d1 */
	if (!parObj.getparfloat("dt", &dt)) dt = ((double) tr.dt)/1000000.0;
	if (!dt) {
		dt = .004;
		warn("dt not set, assumed to be .004");
	}


        parObj.checkpars();


	/* Set up pfa fft */
	nfft = npfaro(nt, LOOKFAC * nt);
	if (nfft >= SU_NFLTS || nfft >= PFA_MAX)
		 throw cseis_geolib::csException("Padded nt=%d--too big", nfft);
	nfby2p1 = nfft/2 + 1;
	d1 = 1.0/(nfft*dt);

	rt = ealloc1float(nfft);
	ct = ealloc1complex(nfby2p1);


	/* Main loop over traces */
	do {
		++ntr;

		/* Load trace into rt (zero-padded) */
		memcpy( (void *) rt, (const void *) tr.data, nt*FSIZE);
		memset( (void *) (rt + nt), 0, (nfft - nt)*FSIZE);

		/* FFT */
		pfarc(1, nfft, rt, ct);

		/* Compute amplitude spectrum */
		tr.data[0] = rcabs(ct[0])/2.0;
		for (i = 1; i < nfby2p1; ++i)  tr.data[i] = rcabs(ct[i]);

		/* Set header values */
		tr.ns = nfby2p1;
		tr.dt = 0;	  /* d1=df is now the relevant step size */
		tr.trid = AMPLITUDE;
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
