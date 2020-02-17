/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUBFILT: $Revision: 1.21 $ ; $Date: 2011/11/12 00:09:00 $	*/

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
std::string sdoc_subfilt =
" 								"
" SUBFILT - apply Butterworth bandpass filter 			"
" 								"
" subfilt <stdin >stdout [optional parameters]			"
" 							        "
" Required parameters:						"
" 	if dt is not set in header, then dt is mandatory	"
" 							        "
" Optional parameters: (nyquist calculated internally)		"
" 	zerophase=1		=0 for minimum phase filter 	"
" 	locut=1			=0 for no low cut filter 	"
" 	hicut=1			=0 for no high cut filter 	"
" 	fstoplo=0.10*(nyq)	freq(Hz) in low cut stop band	"
" 	astoplo=0.05		upper bound on amp at fstoplo 	"
" 	fpasslo=0.15*(nyq)	freq(Hz) in low cut pass band	"
" 	apasslo=0.95		lower bound on amp at fpasslo 	"
" 	fpasshi=0.40*(nyq)	freq(Hz) in high cut pass band	"
" 	apasshi=0.95		lower bound on amp at fpasshi 	"
" 	fstophi=0.55*(nyq)	freq(Hz) in high cut stop band	"
" 	astophi=0.05		upper bound on amp at fstophi 	"
" 	verbose=0		=1 for filter design info 	"
" 	dt = (from header)	time sampling interval (sec)	"
" 							        "
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace subfilt {


/* Credits:
 *	CWP: Dave for bf.c subs and test drivers
 *	CWP: Jack for su wrapper
 *
 * Caveat: zerophase will not do good if trace has a spike near
 *	   the end.  One could make a try at getting the "effective"
 *	   length of the causal filter, but padding the traces seems
 *	   painful in an already expensive algorithm.
 *
 * Trace header fields accessed: ns, dt, trid
 */
/**************** end self doc ***********************************/



segy tr;

void* main_subfilt( void* args )
{
	int zerophase;		/* flag for zero phase filtering	*/
	int locut;		/* flag for low cut filtering		*/
	int hicut;		/* flag for high cut filtering		*/
	float fstoplo;		/* left lower corner frequency		*/
	float fpasslo;		/* left upper corner frequency		*/
	float fpasshi;		/* right lower corner frequency		*/
	float fstophi;		/* right upper corner frequency		*/
	float astoplo;		/* amp at fstoplo			*/
	float apasslo;		/* amp at fpasslo			*/
	float apasshi;		/* amp at fpasshi			*/
	float astophi;		/* amp at fstophi			*/
	int npoleslo;		/* poles in low cut filter		*/
	int npoleshi;		/* poles in high cut filter		*/
	float f3dblo;		/* 3 db point of low cut filter		*/
	float f3dbhi;		/* 3 db point of high cut filter	*/
	float dt;		/* sample spacing			*/
	float nyq;		/* nyquist frequency			*/
	int nt;			/* number of points on input trace	*/
	int verbose;		/* design info flag 			*/
	cwp_Bool seismic;	/* is this seismic data?		*/

	
	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_subfilt );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get info from first trace */ 
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	seismic = ISSEISMIC(tr.trid); 
		 
	if (!seismic)
		warn("input is not seismic data, trid=%d", tr.trid);
	nt = tr.ns;
	if (!parObj.getparfloat("dt", &dt))	dt = ((double) tr.dt)/1000000.0;
	if (!dt) throw cseis_geolib::csException("dt field is zero and not parObj.getparred");
	nyq = 0.5/dt;


	/* Get design frequencies and amplitudes */
	if (!parObj.getparint("verbose", &verbose))	verbose = 0;
	if (!parObj.getparint("zerophase", &zerophase)) zerophase = 1;
	if (!parObj.getparint("locut", &locut))	locut = 1;
	if (!parObj.getparint("hicut", &hicut))	hicut = 1;
	if (!parObj.getparfloat("fstoplo", &fstoplo))	fstoplo = .10 * nyq;
	if (!parObj.getparfloat("fpasslo", &fpasslo))	fpasslo = .15 * nyq;
	if (!parObj.getparfloat("fpasshi", &fpasshi))	fpasshi = .40 * nyq;
	if (!parObj.getparfloat("fstophi", &fstophi))	fstophi = .55 * nyq;
	if (locut) {
		if (fstoplo <= 0.0)      throw cseis_geolib::csException("fstoplo must be positive");
		if (fstoplo > fpasslo)  throw cseis_geolib::csException("fstoplo must be < fpasslo");
	}
	if (hicut) {
		if (fpasshi > fstophi)  throw cseis_geolib::csException("fpasshi must be < fstophi");
		if (fstophi > nyq)  throw cseis_geolib::csException("fstophi must be < nyquist (%f)", nyq);
	}
	if (!parObj.getparfloat("astoplo", &astoplo))	astoplo = .05;
	if (!parObj.getparfloat("apasslo", &apasslo))	apasslo = .95;
	if (!parObj.getparfloat("apasshi", &apasshi))	apasshi = .95;
	if (!parObj.getparfloat("astophi", &astophi))	astophi = .05;
	if (astoplo > apasslo || apasshi < astophi)
		throw cseis_geolib::csException("Bad amplitude parameters");
		
		
	/* Normalize frequencies to [0, 0.5] for bfdesign */
	fstoplo *= dt;
	fpasslo *= dt;
	fstophi *= dt;
	fpasshi *= dt;
	
	
	/* Adapt user frequencies if zerophase selected */
	if (zerophase) {	
		astoplo = sqrt(astoplo);
		apasslo = sqrt(apasslo);
		astophi = sqrt(astophi);
		apasshi = sqrt(apasshi);
	}

	
	/* Use bdesign to make low and high cut filters */
	if (locut) bfdesign(fpasslo,apasslo,fstoplo,astoplo,&npoleslo,&f3dblo);
	if (hicut) bfdesign(fpasshi,apasshi,fstophi,astophi,&npoleshi,&f3dbhi);


	/* Give verbose info if requested */
	if (verbose && locut) {
		if (zerophase) {
			warn("low-cut filter: npoles = %d, 3db point = %f(Hz)",
				2*npoleslo, f3dblo/dt);
		} else {
			warn("low-cut filter: npoles = %d, 3db point = %f(Hz)",
				npoleslo, f3dblo/dt);
		}
	}
	if (verbose && hicut) {
		if (zerophase) {
			warn("high-cut filter: npoles = %d, 3db point = %f(Hz)",
				2*npoleshi, f3dbhi/dt);
		} else {
			warn("high-cut filter: npoles = %d, 3db point = %f(Hz)",
				npoleshi, f3dbhi/dt);

		}
	}

	/* Main loop over traces */
	do {
		/* low-cut (high pass) filter */
		if (locut) {
		    bfhighpass(npoleslo,f3dblo,nt,tr.data,tr.data);
		    if (zerophase) {
			register int i;
		        for (i=0; i<nt/2; ++i) { /* reverse trace in place */
				register float tmp = tr.data[i];
				tr.data[i] = tr.data[nt-1 - i];
				tr.data[nt-1 - i] = tmp;
			}
		        bfhighpass(npoleslo,f3dblo,nt,tr.data,tr.data);
		        for (i=0; i<nt/2; ++i) { /* flip trace back */
				register float tmp = tr.data[i];
				tr.data[i] = tr.data[nt-1 - i];
				tr.data[nt-1 - i] = tmp;
			}
		    }
		}

		/* high-cut (low pass) filter */
		if (hicut) {
		    bflowpass(npoleshi,f3dbhi,nt,tr.data,tr.data);
		    if (zerophase) {
			register int i;
			for (i=0; i<nt/2; ++i) { /* reverse trace */
				register float tmp = tr.data[i];
				tr.data[i] = tr.data[nt-1 - i];
				tr.data[nt-1 - i] = tmp;
			}
			bflowpass(npoleshi,f3dbhi,nt,tr.data,tr.data);
		        for (i=0; i<nt/2; ++i) { /* flip trace back */
				register float tmp = tr.data[i];
				tr.data[i] = tr.data[nt-1 - i];
				tr.data[nt-1 - i] = tmp;
			}
		    }
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
