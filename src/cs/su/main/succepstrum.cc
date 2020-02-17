/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUCCEPSTRUM: $Revision: 1.4 $ ; $Date: 2011/11/16 23:35:04 $	*/

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
#define TWOPI 2.0*PI
#define TINY FLT_EPSILON
#define CTINY cmplx(FLT_EPSILON,FLT_EPSILON)

/*********************** self documentation *****************************/
std::string sdoc_succepstrum =
" SUCCEPSTRUM - Compute the complex CEPSTRUM of a seismic trace 	"
"									"
"  sucepstrum < stdin > stdout					   	"
"									"
" Required parameters:						  	"
"	none								"
" Optional parameters:						  	"
" sign1=1		sign of real to complex transform		"
" sign2=1		sign of complex to real transform		"
" unwrap=1.0	dphase>= PI/unwrap constitutes a wrap in phase		"
"									"
" Notes:								"
" The cepstrum is defined as the fourier transform of the the decibel   "
" spectrum, as though it were a time domain signal.			"
"									"
" CC(t) = FT[ln[T(omega)] ] = FT[ ln|T(omega)| + i phi(omega) ]		"
"	T(omega) = |T(omega)| exp(i phi(omega))				"
"       phi(omega) = unwrapped phase of T(omega)			"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace succepstrum {


/*
 * Author: John Stockwell, Dec 2010
 *   based on sucepstrum.c by:
 * Credits:
 * Balazs Nemeth of Potash Corporation of Saskatchewan Inc. 
 *  given to CWP in 2008
 *
 */
/**************** end self doc ********************************/


#define LOOKFAC 4	/* Look ahead factor for npfaro   */
#define PFA_MAX 720720  /* Largest allowed nfft	   */

/* Segy data constants */
segy tr;				/* SEGY trace */
segy trout;

void rcceps(int sign1, int sign2, float unwrap, int trend, int zeromean, 
			int nt, float *x, float *c);

void* main_succepstrum( void* args )
{
	int nt;		/* number of time samples per trace		*/
	float dt;	/* time sampling interval			*/
	int sign1;	/* sign on real to complex transform		*/
	int sign2;	/* sign on complex to real transform		*/
	float unwrap;	/* unwrapping divisor 				*/

	int trend;	/* =1 remove trend				*/
	int zeromean;	/* =0 assume phase(0) =0.0, =1 assume zero mean */
	
	/* hook up getpars */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_succepstrum );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	
	/* get parameters */
	if (!parObj.getparint("sign1",&sign1)) 	sign1=1;
	if (!parObj.getparint("sign2",&sign2)) 	sign2=1;
	if (!parObj.getparfloat("unwrap",&unwrap)) 	unwrap=1.0;
	if (!parObj.getparint("trend",&trend)) 	trend=1;
	if (!parObj.getparint("zeromean",&zeromean)) 	zeromean=0;
	

        parObj.checkpars();

	/* get information from the first header */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;
	dt = (double)tr.dt/1000000.0;
		
	/* loop over traces */
	do {

		/* compute the cepstrum or compute	*/
		rcceps(sign1,sign2,unwrap,trend,zeromean,nt,tr.data,tr.data);	

       		/* Set header values--npfaro makes nfft even */
		tr.trid = FUNPACKNYQ;

			
		su2cs->putTrace(&tr); 
	} while(cs2su->getTrace(&tr));
	
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

void rcceps(int sign1,int sign2, float unwrap, int trend, int zeromean, int nt, 
		float *x,float *c)
/***********************************************************************
rcceps - compute the complex cepstrum  of a signal
************************************************************************
Input:
nt	number of time samples
x	input data
Output:
c	output data complex
************************************************************************
Notes:
The complex cepstrum is defined as 
C(t) = FT[ log[ F(omega) ] ]

F(omega) = FT[ F(t) ] 

F(omega) is a complex variable that can be written in polar form as:
F(omega) = | F(omega) | exp( i phi(omega) )

where phi(omega) is the unwrapped phase of F(omega)

log[F(omega)] = log|F(omega)| + i phi(omega)

C(t) = FT[  log|F(omega)| + i phi(omega) ]

************************************************************************
Author: CWP, John Stockwell, Dec 2010 based on rceps by:
Credits: Balasz Nemeth, Potash Corporation, Saskatchewan  c. 2008
***********************************************************************/
{	
	int nfftc;		/* padded size of nfft		*/
	int nf;			/* number of frequencies	*/
	float snfftc;		/* scaling factor 		*/
	complex *w=NULL;	/* frequency domain transform	*/
	float *a=NULL;		/* amplitude			*/
	float *p=NULL;		/* phase 			*/
	int iw;			/* loop counter for frequency	*/
	
	
	/* Set up pfa fft */
	nfftc = npfao(nt,LOOKFAC*nt); 
	if (nfftc >= SU_NFLTS || nfftc >= PFA_MAX)
			throw cseis_geolib::csException("Padded nt=%d--too big", nfftc);
	nf = nfftc/2 + 1;
	snfftc=1.0/nfftc;

	/* allocate space */
	w = ealloc1complex(nf);
	a = ealloc1float(nfftc);
	p = ealloc1float(nfftc);
		
	memset( (void *) &a[nt], 0, (nfftc-nt)*FSIZE);
	memcpy( (void *) a, (const void *) x, nt*FSIZE);
		
	/* FFT */			
	sscal(nt,snfftc,a,1);
	pfarc(sign1, nfftc,a,w);
	

	/* find amplitudes and phases */
	for(iw=0;iw<nf;++iw) {
		a[iw] = rcabs(w[iw]);
		p[iw] = atan2(w[iw].i,w[iw].r);
	}

	/* unwrap the phase */
	simple_unwrap_phase(nf, trend, zeromean, unwrap,p);
	
	/* take the log of the amplitude */
	for(iw=0;iw<nf;++iw) {
		if(!CLOSETO(a[iw],0.0)) {
			w[iw].r = (float)log((double)a[iw]);
			w[iw].i = p[iw];
		} else {
			w[iw].r=0.0;
			w[iw].i=0.0;
		}
		
	}
	pfacc(sign2, nf,w);

	for(iw=0;iw<nf;++iw) {
		c[2*iw]   = w[iw].r;
		c[2*iw+1] = w[iw].i;
	}

	
	free1float(a);
	free1float(p);
	free1complex(w);

}


} // END namespace
