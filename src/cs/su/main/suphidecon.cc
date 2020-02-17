/* SUPHIDECON: $Revision: 1.8 $ ; $Date: 2011/12/21 23:23:59 $	*/


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

/*********************** self documentation *****************************/
std::string sdoc_suphidecon =
"									"
" SUPHIDECON - PHase Inversion Deconvolution				"
"									"
"    suphidecon < stdin > stdout					"
"									"
" Required parameters:						  	"
"	none							   	"
" Optional parameters:							"
" ... time range used for wavelet extraction:			   	"
" tm=-0.1	Pre zero time (maximum phase component )		"
" tp=+0.4	Post zero time (minimum phase component + multiples)    "
" percpad=50	percentage padding for nt prior to cepstrum calculation	"
"									"
" pnoise=0.001	Pre-withening (assumed noise to prevent division by zero)"
"									"
" Notes:								"
" The wavelet is separated from the reflectivity and noise based on	"
" their different 'smoothness' in the pseudo cepstrum domain.		"
" The extracted wavelet also includes multiples. 			"
" The wavelet is reconstructed in frequency domain, end removed		" 
" from the trace. (Method by Lichman and Northwood, 1996.)		"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suphidecon {


/*
 *
 * Credits: Potash Corporation, Saskatechwan  Balasz Nemeth 
 * given to CWP by Potash Corporation 2008 (originally as supid.c)
 *
 * Reference:
 * Lichman,and Northwood, 1996; Phase Inversion deconvolution for
 * long and short period multiples attenuation, in
 * SEG Deconvolution 2, Geophysics reprint Series No. 17
 * p. 701-718, originally presented at the 54th EAGE meeting, Paris,
 * June 1992, revised March 1993, revision accepted September 1994.
 * 
 *
 */

/**************** end self doc ********************************/

#define LOOKFAC 4	 /* Look ahead factor for npfaro   */
#define PFA_MAX 720720  /* Largest allowed nfft	   */

/* Segy data constants */
segy tr;				/* SEGY trace */

/* function prototype of routine used internally */
int computePseudoCepstrum(int *nt,float *percpad, float *x,complex *c,int init);

void* main_suphidecon( void* args )
{
	int nt;			/* number of time samples per trace	*/
	int it;			/* counter				*/
	float dt;		/* time sampling interval		*/
	int nfft;		/* number of points in fft		*/
	float tm;		/* maximum phase component		*/
	float tp;		/* minimum phase component + multiples	*/
	int itmin,itmax;	/* min and max time samples 		*/
	int nwind;		/* number of windows 			*/
	float *window=NULL;	/* hanning window function		*/
	int index;		/* index				*/
	
	float pnoise;		/* prewitening to prevent division by 0	*/
	float percpad;		/* pad nt by percpad percentage prior	*/
				/* to npfa padding calculation		*/
	
	complex *c=NULL;	/* complex array			*/
	complex *ct=NULL;	/* complex array			*/
	float *rt=NULL;		/* input (real) data			*/

	
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suphidecon );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	
	if (!parObj.getparfloat("tm", &tm)) tm=-0.1;
	if (!parObj.getparfloat("tp", &tp)) tp=0.4;
	if (!parObj.getparfloat("pnoise",&pnoise)) pnoise=0.1;
	
	/* get information from the first header */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;
	dt = (double)tr.dt/1000000.0;
	
	itmin=-tp/dt;
	itmax=-tm/dt;

	nwind = itmax-itmin+1;

	/* allocate space for hanning window */
	window = ealloc1float(nwind);

	/* compute hanning window */
	hanningnWindow(nwind,window);
	
	/* get nfft as required for cepstrem */
	nfft = computePseudoCepstrum(&nt,&percpad,rt,c,1);	

	/* allocate space */
	c = ealloc1complex(nfft);
	ct = ealloc1complex(nfft);
	rt = ealloc1float(nfft);

		
	do {

			memset((void *)rt, 0, nfft*FSIZE);
			memset((void *)ct, 0, nfft*sizeof(complex));
			memset((void *)c, 0, nfft*sizeof(complex));

                       /* copy traces from data */
                        memcpy(( void *) rt, (const void *) tr.data, nt*FSIZE);

			computePseudoCepstrum(&nt,&percpad,rt,c,0);	
			
			/* non-symetrical linear filter 
			   from itmin to itmax; Hanning window */
			
		/* wrap array so memory is aligned from itmin to itmax */
			wrapArray(c,nfft,sizeof(complex),itmax+1);
			
			/* Apply hanning window */
			for(it=0,index=0;it<nwind;++it,++index)
				c[it] = crmul(c[it],window[index]);
				
			/* Unwrapp array */
			wrapArray(c,nfft,sizeof(complex),-(itmax+1));
			

			pfacc(-1, nfft,c); 
			sscal(nfft*2,1.0/nfft,(float*)c,1); 
			
			
			/* remove log ; exponentiate */
			{ complex tmpc;
				for(it=0;it<nfft;it++) {
					tmpc = cwp_cexp(c[it]);
					if(isfinite(tmpc.r) && isfinite(tmpc.i)) {
						c[it] = tmpc;
					} else {
						c[it] = cmplx(0.0,0.0);
					}	
				}	
			}
			
			
			/* Shift it to the central location */
		/*	for(it=1;it<nfft;it+=2) {
				c[it].r = -c[it].r;
				c[it].i = -c[it].i;
			  }
		*/


			
			/* Deconvolution */
			/* Fourier transform trace */
			pfarc(-1, nfft,rt,ct);
			sscal(nfft*2,1.0/nfft,(float*)ct,1); 
			
			/* Pre-whiten */
			{ float Spsum=0.0;
			  int nf=nfft/2+1;
			
				for(it=0;it<nf;it++)
					Spsum += rcabs(c[it]);
				pnoise = pnoise/10000.0*(Spsum/nf);
			
				/* create full spectrum */
				for(it=0;it<nf-2;it++) {
					ct[nf+it].r = ct[nf-it-2].r;
					ct[nf+it].i = -ct[nf-it-2].i;
				}
			}
			
			/* Deconvolve */
			for(it=0;it<nfft;it++)
				ct[it] = cdiv(ct[it],cmplx(c[it].r+pnoise,c[it].i));
			pfacr(1, nfft,ct,rt);
			
			for(it=0;it<nfft;it++) {
				tr.data[it] = rt[it];
			}
			
			su2cs->putTrace(&tr);
			
	} while(cs2su->getTrace(&tr));
	
/*
	free1complex(c);
	free1complex(ct);
	free1float(rt);
	
*/
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

int computePseudoCepstrum(int *nt,float *percpad,float *x,complex *c,int init)
/*********************************************************************
 computePseudoCepstrum  - compute the pseudo cepstrum
*/
{	
	int nfftc;		/* nfft for the cepstrum */
	int nf;			/* number of frequencies */
	float snfftc;		/* 1.0/nfftc 		*/		
	complex *w=NULL;
	float *a=NULL;
	float *p=NULL;
	int iw;
	int ntp;
	
	
	/* Set up pfa fft */
	ntp = NINT(*nt*(1.0+ *percpad/100.0));
	nfftc = npfao(ntp,LOOKFAC*ntp); 
	if (nfftc >= SU_NFLTS || nfftc >= PFA_MAX)
		 	throw cseis_geolib::csException("Padded nt=%d--too big", nfftc);
	nf = nfftc/2 + 1;
	snfftc=1.0/nfftc;
	
	if (init)
		return(nfftc);
		
	w = ealloc1complex(nfftc);
	a = ealloc1float(nfftc);
	p = ealloc1float(nfftc);
		
	memset( (void *) &a[*nt],0, (nfftc-*nt)*FSIZE);
	memcpy( (void *) a, (const void *) x, *nt*FSIZE);
		
	/* FFT */			
	sscal(*nt,snfftc,a,1);
	pfarc(-1, nfftc,a,w);
	
	for(iw=0;iw<nf;iw++) {
		a[iw] = w[iw].r;
	}
	
	hilbert(nf,a,p);

	for(iw=0;iw<nf;iw++) {
		w[iw].i=(-1.0/PI)*p[iw];
	}
	
	for(iw=0;iw<nf;iw++) {
		if(!CLOSETO(rcabs(w[iw]),0.0)) {
			/*w[iw] = crmul(cwp_clog(w[iw]),0.5); */
			w[iw] = cwp_clog(w[iw]);
		} else {
			w[iw].r=0.0;
			w[iw].i=0.0;
		}
	}
	
	/* prepare it for in place complex transform */

	for(iw=0;iw<nf-2;iw++) {
		w[nf+iw].r = w[nf-iw-2].r;
		w[nf+iw].i = -w[nf-iw-2].i;
	}
	
	/* Shift it to the central location */
	/*for(iw=0;iw<nfftc;iw+=2) {
		w[iw].r = -w[iw].r;
		w[iw].i = -w[iw].i;
	} */


	pfacc(1, nfftc,w);
					
	for(iw=0;iw<nfftc;iw++) {
		c[iw] = w[iw];
	/*	c[iw] = cmul(w[iw],cmplx(-1.0,-1.0)); */
	}
	
/*
	free1float(a);
	free1float(p);
	free1complex(w);
*/
	
	return(nfftc);
}


} // END namespace
