/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SULOG: $Revision: 1.19 $ ; $Date: 2011/11/16 23:21:55 $		*/

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
std::string sdoc_sulog =
"								"
" SULOG -- time axis log-stretch of seismic traces		"
"								"
" sulog [optional parameters] <stdin >stdout 			"
"								"
" Required parameters:						"
"	none				 			"
"								"
" Optional parameters:						"
"	ntmin= .1*nt		minimum time sample of interest	"
"	outpar=/dev/tty		output parameter file, contains:"
"				number of samples (nt=)		"
"				minimum time sample (ntmin=)	"
"				output number of samples (ntau=)"
"	m=3			length of stretched data	"
"				is set according to		"
"					ntau = nextpow(m*nt)	"
"	ntau= pow of 2		override for length of stretched"
"				data (useful for padding zeros	"
"				to avoid aliasing)		"
"								"
" NOTES:							"
"	ntmin is required to avoid taking log of zero and to 	"
"	keep number of outsamples (ntau) from becoming enormous."
"        Data above ntmin is zeroed out.			"
"								"
"	The output parameters will be needed by suilog to 	"
"	reconstruct the original data. 				"
"								"
" EXAMPLE PROCESSING SEQUENCE:					"
"		sulog outpar=logpar <data1 >data2		"
"		suilog par=logpar <data2 >data3			"
"								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sulog {


/* Credits:
 *	CWP: Shuki, Chris
 *
 * Caveats:
 * 	Amplitudes are not well preserved.
 *
 * Trace header fields accessed: ns, dt
 * Trace header fields modified: ns, dt
 */
/**************** end self doc ***********************************/


/* prototypes of functions used internally */
int nextpower(int p, int n); /* function for padding to power of 2 */
void stretch(float *q, float *p, float *w, int *it, int lq, int nw);
void lintrp(float *q, float *w, int *it, int lp, int lq);

segy tr;

void* main_sulog( void* args )
{
	float *buf;	/* temporary repository of log stretched data 	*/
	float dt;	/* input time sampling interval			*/
	float dtau;	/* output time sampling interval 		*/
	float *t;	/* fractional sample number on input data	*/
	float *w;	/* Interpolation weights 			*/
	int *it;	/* Interpolation locations			*/
	int itau;	/* tau sample counter				*/
	int m;		/* m*nt=ntau 					*/
	int nt;		/* time samples on input data			*/
	int ntau;	/* pow of 2 >= ntautemp  ==  samples in outdata	*/
	size_t ntausize;/*  ... in bytes				*/
	int ntautemp;	/* tau sample corresponding to nt	 	*/
	int ntmin;	/* minimum input time of interest		*/
	int nw;		/* Number of interpolation weights (2=linear)	*/
	cwp_String outpar;/* name of file holding output parfile	*/
	FILE *outparfp;	/* ... its file pointer				*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sulog );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get information from the first header */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;
	dt = ((double) tr.dt)/1000000.0;

	/* Minimum sample of interest (data above tmin is lost) */
	if (!parObj.getparint("ntmin", &ntmin))	ntmin = 0.1 * nt;
	if (!parObj.getparstring("outpar", &outpar))	outpar = "/dev/tty";

	/* Open file to save parameters */
	outparfp = efopen(outpar, "w");

	/* 2 weights for linear interpolation */
	nw = 2;
	
	/* fix length of stretched data (but will be padded) */
	if (!parObj.getparint("m", &m))	 m = 3;
	ntautemp = m*nt;  

	/* find dtau, dt in log argument cancels */ 
	dtau = log((float) nt/ntmin) / ntautemp;
	if (!parObj.getparint("ntau", &ntau))	 ntau = nextpower(2, ntautemp);
        parObj.checkpars();
	CHECK_NT("ntau",ntau);
	ntausize = ntau * FSIZE;
	
	/* Put out enough parameters to reconstruct original data */
	/* (ntau for reassurance to user, also in tr.ns of output) */
	(void) fprintf(outparfp,
		       "nt=%d ntmin=%d dt=%g ntau=%d\nt=%d n1=%d d1=%g", 
		       nt, ntmin, dt*1000000.0, ntau, nt, nt, dt*1000000.0);

	/* Allocate space for stretching operation */
	t   = ealloc1float(ntau);
	it  = ealloc1int(ntau);
	w   = ealloc1float(nw * ntau);
	buf = ealloc1float(ntau);


/* 	The log-stretch from 't' to 'tau' is given mathematically by  	
 *  								
 * 		tau = log( t / tmin ) + taumin           
 *  							
 * 	taumin is arbitrary and taken to be taumin=0	
 */
  					
	/* Calculate fractional t-sample that each tau-sample maps to;  */
	/* [ itau=0 --> (float)itmin ..&.. itau=ntau --> (float)itmax ] */
	for (itau = 0; itau < ntau; itau++) {
		t[itau] = (float) ntmin * exp((float) itau*dtau);
	}

	/* Calculate the linear interpolation coefficients */
	lintrp (t, w, it, nt, ntau);

	/* Main loop over traces */
	do {
		/* Perform the stretch; put new data into buf */
		stretch (buf, tr.data, w, it, ntau, nw);

		/* Overwrite the segy data */
		memcpy( (void *) tr.data, (const void *) buf, ntausize);

		tr.ns = ntau;
		tr.dt = dtau*1000000.;

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

void stretch(float *q, float *p, float *w, int *it, int lq, int nw)
/*
 *  General coordinate stretch with predetermined coefficients
 *
 *         NW-1
 * Q(T) =  SUM W(T,J)*P(IT(T)), FOR T=0,LQ-1
 *         J=0
 */
{
	int j, i;

	for (i = 0; i < lq; i++) {
		q[i] = 0.0;
		for (j=0; j<nw; j++) {
			q[i] += w[i*nw+j] * p[it[i]+j];
		}
	}
	return;
}

void lintrp (float *q, float *w, int *it, int lp, int lq)
{
	int i;
	float delta;

	for (i = 0; i < lq; i++) {
		if (q[i] >= 0.0 && q[i] < lp - 1) {
			it[i] = q[i]; 
			delta = q[i] - it[i];
			w[i*2] = 1.0 - delta;
			w[i*2+1] = delta;
		} else {
			it[i] = 0;
			w[i*2] = 0.0;
			w[i*2+1] = 0.0;
		}
	}
	return;
}

int nextpower(int p, int n)
{
	int nn;
	if (!n) return 0;
	for (nn = 1; nn < n; nn *= p);
	return nn;
}

} // END namespace
