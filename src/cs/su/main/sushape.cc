/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSHAPE: $Revision: 1.15 $ ; $Date: 2011/11/16 17:47:47 $		*/

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

/*********************** self documentation ******************************/
std::string sdoc_sushape =
" 									"
" SUSHAPE - Wiener shaping filter					"
" 									"
"  sushape <stdin >stdout  [optional parameters]			"
" 									"
" Required parameters:							"
" w=		vector of input wavelet to be shaped or ...		"
" ...or ... 								"
" wfile=        ... file containing input wavelet in SU (SEGY trace) format"
" d=		vector of desired output wavelet or ...			"
" ...or ... 								"
" dfile=        ... file containing desired output wavelet in SU format	"
" dt=tr.dt		if tr.dt is not set in header, then dt is mandatory"
" 									"
" Optional parameters:							"
" nshape=trace		length of shaping filter			"
" pnoise=0.001		relative additive noise level			"
" showshaper=0		=1 to show shaping filter 			"
" 									"
" verbose=0		silent; =1 chatty				"
" 									"
"Notes:									"
" 									"
" Example of commandline input wavelets: 				"
"sushape < indata  w=0,-.1,.1,... d=0,-.1,1,.1,... > shaped_data	"
" 									"
"sushape < indata  wfile=inputwavelet.su dfile=desire.su > shaped_data	"
" 									"
" To get the shaping filters into an ascii file:			"
" ... | sushape ... showwshaper=1 2>file | ...   (sh or ksh)		"
" (... | sushape ... showshaper=1 | ...) >&file  (csh)			"
" 									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sushape {


/* Credits:
 *	CWP: Jack Cohen
 *	CWP: John Stockwell, added wfile and dfile  options
 *
 * Trace header fields accessed: ns, dt
 * Trace header fields modified: none
 *
 */
/**************** end self doc *******************************************/


#define PNOISE	0.001


segy intrace, outtrace;
segy dtr, wtr;

void* main_sushape( void* args )
{
	int nt;			/* number of points on trace		*/

	float dt;		/* time sample interval (sec)		*/
	float *shaper;		/* shaping filter coefficients		*/
	float *spiker;		/* spiking decon filter (not used)	*/
	float *w;		/* input wavelet			*/

	int nw;			/* length of input wavelet in samples	*/
	float *d;		/* desired output wavelet		*/

	int nd;			/* length of desired wavelet in samples	*/
	int nshape;		/* length of shaping filter in samples	*/

	float pnoise;		/* pef additive noise level		*/
	float *crosscorr;	/* right hand side of Wiener eqs	*/
	float *autocorr;	/* vector of autocorrelations		*/
	int showshaper;		/* flag to display shaping filter	*/
        float f_zero=0.0;       /* zero valued item for comparison      */

	cwp_String wfile="";	/* input wavelet file name		*/
	cwp_String dfile="";	/* desired output wavelet file name	*/
	FILE *wfp;		/* input wavelet file pointer 		*/
	FILE *dfp;		/* desired wavelet file pointer		*/
	int verbose=0;		/* =0 silent; =1 chatty			*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sushape );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get info from first trace */ 
	if (!cs2su->getTrace(&intrace)) throw cseis_geolib::csException("can't get first trace");
	nt = intrace.ns;
	dt = intrace.dt/1000000.0;	if (!dt) CSMUSTGETPARFLOAT ("dt", &dt);


	/* Get parameters */
	if (!parObj.getparint("showshaper",  &showshaper))	showshaper = 0;
	if (!parObj.getparint("nshape",  &nshape))		nshape = nt;
	if (!parObj.getparfloat("pnoise",  &pnoise))		pnoise = PNOISE;
	if (!parObj.getparint("verbose", &verbose))		verbose = 0;

	/* Open dfile and wfile if they have been parObj.getparred */
	parObj.getparstring("dfile",&dfile);	
	parObj.getparstring("wfile",&wfile);	

	if ((*dfile=='\0')) { /* if no dfile, then get from command line */
		if (!(nd = parObj.countparval("d")))
			throw cseis_geolib::csException("must specify d= desired wavelet");
		d = ealloc1float(nd);	parObj.getparfloat("d", d);

	} else { /* read from dfile  */

                if((dfp=fopen(dfile,"r"))==NULL)
                        throw cseis_geolib::csException("cannot open dfile=%s\n",dfile);

        	if (!fgettr(dfp,&dtr))  throw cseis_geolib::csException("can't get input wavelet");
        		nd = (int) dtr.ns;
		d = ealloc1float(nd);
		memcpy((void *) d, (const void *) dtr.data, nd*FSIZE);
	}
		
	if ((*wfile=='\0')) { /* then get w from command line */
		if (!(nw = parObj.countparval("w")))
			throw cseis_geolib::csException("must specify w= desired wavelet");
		w = ealloc1float(nw);	parObj.getparfloat("w", w);

	} else { /* read from wfile  */

                if((wfp=fopen(wfile,"r"))==NULL)
                        throw cseis_geolib::csException("cannot open wfile=%s\n",wfile);

        	if (!fgettr(wfp,&wtr))  throw cseis_geolib::csException("can't get desired output wavelet");
        		nw = (int) wtr.ns;
		w = ealloc1float(nw);
		memcpy((void *) w, (const void *) wtr.data, nw*FSIZE);
	}

        parObj.checkpars();

	/* Get shaping filter by Wiener-Levinson */
	shaper	  = ealloc1float(nshape);
	spiker 	  = ealloc1float(nshape);	/* not used */
	crosscorr = ealloc1float(nshape);
	autocorr  = ealloc1float(nshape);
	xcor(nw, 0, w, nw, 0, w, nshape, 0, autocorr);  /* for matrix */
	xcor(nw, 0, w, nd, 0, d, nshape, 0, crosscorr); /* right hand side */
        if (CLOSETO(autocorr[0],f_zero))  throw cseis_geolib::csException("can't shape with zero wavelet");
	autocorr[0] *= (1.0 + pnoise);			/* whiten */
	stoepf(nshape, autocorr, crosscorr, shaper, spiker);
		

	/* Show shaper on request */
	if (showshaper) {
		register int i;
		if (verbose) warn("Shaping filter:");
		for (i = 0; i < nshape; ++i)
			fprintf(stderr, "%10g%c", shaper[i],
				(i%6==5 || i==nshape-1) ? '\n' : ' ');
	}



	/* Main loop over traces */
	do {
		/* Center and convolve shaping filter with trace */
		convolve_cwp(nshape, (nw-nd)/2, shaper,
		     nt, 0, intrace.data, 
                     nt, 0, outtrace.data);        


		/* Output filtered trace */
		memcpy( (void *) &outtrace, (const void *) &intrace, HDRBYTES);
		su2cs->putTrace(&outtrace);

	} while (cs2su->getTrace(&intrace));


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
