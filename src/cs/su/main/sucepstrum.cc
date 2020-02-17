/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUCEPSTRUM: $Revision: 1.7 $ ; $Date: 2011/11/16 23:35:04 $	*/

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
std::string sdoc_sucepstrum =
" SUCEPSTRUM - Compute the CEPSTRUM of a seismic trace or compute the	"
"		minimum phase reconstruction of a trace			"
"  sucepstrum < stdin > stdout					   	"
"									"
" Required parameters:						  	"
"	none								"
" Optional parameters:						  	"
" sign1=1		sign of real to complex transform		"
" sign2=1		sign of complex to real transform		"
"   mpr=0		=1 minimum phase reconstructed version of	" 
"			    original signal				"
"   centered=0		=1 symmetrical output			 	"
" Notes:								"
" The cepstrum is defined as the fourier transform of the the decibel   "
" spectrum, as though it were a time domain signal.			"
" C(t) = FFT(ln(|T(w)|)}						"
"	|T(w)| = amplitude spectrum of the trace			"
" For mpr=1, the input traces are converted to the minimum phase form.	"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sucepstrum {


/*
 * Author: Balazs Nemeth of Potash Corporation of Saskatchewan Inc. 
 *  given to CWP in 2008
 *
 */
/**************** end self doc ********************************/

#define LOOKFAC 4	/* Look ahead factor for npfaro   */
#define PFA_MAX 720720  /* Largest allowed nfft	   */

/* Segy data constants */
segy tr;				/* SEGY trace */
segy trout;

void rceps(int sign1, int sign2, int nt, int mph,float *x,float *c);

void* main_sucepstrum( void* args )
{
	int nt;		/* number of time samples per trace		*/
	float dt;	/* time sampling interval			*/
	int centered;	/* flag; =1 symmetric output 			*/
	int mpr;	/* flag: =1 perform minimum phase reconstruction*/
	int sign1;	/* sign on real to complex transform		*/
	int sign2;	/* sign on complex to real transform		*/
	
	/* hook up getpars */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sucepstrum );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	
	/* get parameters */
	if (!parObj.getparint("centered", &centered)) centered=0;
	if (!parObj.getparint("mpr", &mpr))		mpr=0;
	if (!parObj.getparint("sign1",&sign1)) 	sign1=1;
	if (!parObj.getparint("sign2",&sign2)) 	sign2=1;
	

        parObj.checkpars();

	/* get information from the first header */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;
	dt = (double)tr.dt/1000000.0;
		
	/* loop over traces */
	do {

			/* compute the cepstrum or compute	*/
			/* minimum phase reconstruction		*/
			rceps(sign1,sign2,nt,mpr,tr.data,tr.data);	
			
			/* computation if centered output desired */
			if(centered) {
				{ int n,it;
					if(!ISODD(nt)) {
						n=nt+1;
					} else {
						n=nt;
					}
					memcpy((void *) &trout,(const void *) &tr,HDRBYTES);
					memcpy((void *) &trout.data[n/2],(const void *) tr.data,n/2*FSIZE);
					for(it=0;it<n/2;it++) 
						trout.data[it] = tr.data[n/2-it];
					trout.ns=n;
					trout.f1=-n/2*dt;
					trout.d1=dt;
					su2cs->putTrace(&trout);
				}
			} else {
				su2cs->putTrace(&tr);
			}
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

void rceps(int sign1,int sign2,int nt, int mph,float *x,float *c)
/***********************************************************************
rceps - compute the cepstrum or the minimum phase reconstruction of a signal
        using the cepstrum.
************************************************************************
Input:
nt	number of time samples
mph	minimum phase reconstruction flag
x	input data
Output:
c	output data
************************************************************************
Author: Balasz Nemeth, Potash Corporation, Saskatchewan  c. 2008
***********************************************************************/
{	
	int nfftc;
	int nf;
	float snfftc;
	complex *w=NULL;
	float *a=NULL;
	float *p=NULL;
	int iw;
	int ntp;
	
	
	/* Set up pfa fft */
	ntp = NINT(nt);
	nfftc = npfao(ntp,LOOKFAC*ntp); 
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
	

	for(iw=0;iw<nf;iw++) {
		a[iw] = rcabs(w[iw]);
		p[iw] = atan2(w[iw].i,w[iw].r);
	}
	
	for(iw=0;iw<nf;iw++) {
		if(!CLOSETO(a[iw],0.0)) {
			w[iw].r = (float)log((double)a[iw]);
			w[iw].i = 0.0;
		} else {
			w[iw].r=0.0;
			w[iw].i=0.0;
		}
	}
	pfacr(sign2, nfftc,w,c);

	if(mph) {

		a[0] = c[0];
		if(nt%2) {
			for(iw=1;iw<nt/2;iw++)
				a[iw] = 2.0*c[iw];
			a[nt/2] = c[nt/2];
			memset( (void *) &a[nt/2+1],0,(nt/2-1)*FSIZE);
		} else {
			for(iw=1;iw<nt/2;iw++)
				a[iw] = 2.0*c[iw];
			memset( (void *) &a[nt/2],0,(nt/2-1)*FSIZE);
		}
		
		memset( (void *) &a[nt], 0, (nfftc-nt)*FSIZE);

		sscal(nt,snfftc,a,1);
		pfarc(sign2, nfftc,a,w);
		
		for(iw=0;iw<nf;iw++)
			if(!CLOSETO(rcabs(w[iw]),0.0)) {
				w[iw] = cwp_cexp(w[iw]);
			} else {
				w[iw].r=0.0;
				w[iw].i=0.0;
			}
		
		pfacr(sign2, nfftc,w,c);
	}
	
	free1float(a);
	free1float(p);
	free1complex(w);
}


} // END namespace
