/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUPLANE: $Revision: 1.20 $ ; $Date: 2011/11/12 00:40:42 $	*/

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
std::string sdoc_suplane =
"								"
"SUPLANE - create common offset data file with up to 3 planes	"
"								"
"suplane [optional parameters] >stdout	 			"
"								"
"Optional Parameters:						"
" npl=3			number of planes			"
" nt=64 		number of time samples			"
" ntr=32		number of traces			"
" taper=0		no end-of-plane taper			"
"			= 1 taper planes to zero at the end	"
" offset=400 		offset					"
" dt=0.004	 	time sample interval in seconds		"
"...plane 1 ...							"
"	dip1=0		dip of plane #1 (ms/trace)		"
" 	len1= 3*ntr/4	HORIZONTAL extent of plane (traces)	"
"	ct1= nt/2	time sample for center pivot	 	"
"	cx1= ntr/2	trace for center pivot			"
"...plane 2 ...							"
"	dip2=4		dip of plane #2 (ms/trace)		"
"	len2= 3*ntr/4	HORIZONTAL extent of plane (traces)	"
"	ct2= nt/2	time sample for center pivot 		"
"	cx2= ntr/2	trace for center pivot			"
"...plane 3 ...							"
"	dip3=8		dip of plane #3 (ms/trace)		"
"	len3= 3*ntr/4	HORIZONTAL extent of plane (traces)	"
"	ct3= nt/2	time sample for center pivot		"
"	cx3= ntr/2	trace for center pivot			"
"								"
" liner=0	use parameters					"
"			= 1 parameters set for 64x64 data set   "
"			with separated dipping planes.		"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suplane {


/* Credits:
 *	CWP: Chris Liner
 *
 * Trace header fields set: ns, dt, offset, tracl
 */
/**************** end self doc ***********************************/


#define NT	64
#define NTR	32
#define DT	0.004
#define OFF	400
#define NPL	3


segy tr;

void* main_suplane( void* args )
{
	float dip1;		/* time-dip of plane 1 (ms/trace)	*/
	float dip2;		/* time-dip of plane 2 (ms/trace)	*/
	float dip3;		/* time-dip of plane 3 (ms/trace)	*/
	float dt;		/* time sample interval in seconds 	*/
	float eps;		/* fit - itless (for linear interpol)	*/
	float fit;		/* fractional time sample		*/
	float offset;		/* constant offset for header 		*/
	int ct1,cx1;		/* center of plane 1 (sample and trace)	*/
	int ct2,cx2;		/* center of plane 2 (sample and trace)	*/
	int ct3,cx3;		/* center of plane 3 (sample and trace)	*/
	int itless;		/* sample above plane intersection	*/
	int itmore;		/* sample below plane intersection 	*/
	int itr;		/* trace counter			*/
	int len1;		/* HORIZ extent of plane 1		*/
	int len2;		/* HORIZ extent of plane 2		*/
	int len3;		/* HORIZ extent of plane 3		*/
	int liner;		/* flag for special output section	*/
	float msps;		/* milliseconds per sample 		*/
	int npl;		/* number of planes 			*/
	int nt;			/* time samples in outdata		*/
	int ntr;		/* traces in outdata			*/
	int tfe1;		/* traces-from-end of plane 1 (for taper) */
	int tfe2;		/* traces-from-end of plane 2 (for taper) */
	int tfe3;		/* traces-from-end of plane 3 (for taper) */
	int taper;		/* flag to taper plane ends to zero	*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suplane );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* set parameters and fill header fields */
	nt = NT;	parObj.getparint("nt", &nt);
	CHECK_NT("nt",nt);				tr.ns = nt;
	ntr = NTR;	parObj.getparint("ntr", &ntr);
	dt = DT;	parObj.getparfloat("dt", &dt);		tr.dt = NINT(dt*1000000.);
	offset=OFF;	parObj.getparfloat("offset", &offset);	tr.offset = offset;
	npl = NPL;	parObj.getparint("npl", &npl);
	taper = 0;	parObj.getparint("taper", &taper);

	/* set defaults and/or get parameters for plane 1 */
	dip1 = 0;	parObj.getparfloat("dip1", &dip1);
	len1 = 3*ntr/4;	parObj.getparint("len1", &len1);
	ct1 = nt/2;	parObj.getparint("ct1", &ct1);		ct1 -= 1;
	cx1 = ntr/2;	parObj.getparint("cx1", &cx1);		cx1 -= 1;

	/* set defaults and/or get parameters for plane 2 */
	dip2 = 4;	parObj.getparfloat("dip2", &dip2);
	len2 = 3*ntr/4;	parObj.getparint("len2", &len2);
	ct2 = nt/2;	parObj.getparint("ct2", &ct2);		ct2 -= 1;
	cx2 = ntr/2;	parObj.getparint("cx2", &cx2);		cx2 -= 1;	

	/* set defaults and/or get parameters for plane 3 */
	dip3 = 8;	parObj.getparfloat("dip3", &dip3);
	len3 = 3*ntr/4;	parObj.getparint("len3", &len3);
	ct3 = nt/2;	parObj.getparint("ct3", &ct3);		ct3 -= 1;
	cx3 = ntr/2;	parObj.getparint("cx3", &cx3);		cx3 -= 1;	

	/* check if user wants the special output specified */
        /* by liner=1; if so, set parameters accordingly    */
	liner = 0;	parObj.getparint("liner", &liner);
	if (liner == 1) {
		nt = 64;	tr.ns = nt;
		ntr = 64;
		npl = 3;	

		dip1 = 0;
		len1 = ntr/4;	
		ct1 = nt/2;	ct1 -= 1;		
		cx1 = 3*ntr/4;	cx1 -= 1;

		dip2 = 4;
		len2 = ntr/4;	
		ct2 = nt/2;	ct2 -= 1;
		cx2 = ntr/2;	cx2 -= 1;	

		dip3 = 8;
		len3 = ntr/4;	
		ct3 = nt/2;	ct3 -= 1;
		cx3 = ntr/4;	cx3 -= 1;	
	}

	/* calculate milliseconds per sample */
	msps = dt*1000.0;	

	tfe1 = 0; tfe2 = 0; tfe3 = 0;
	for (itr = 0; itr < ntr; itr++) {
		memset( (void *) tr.data, 0, nt * FSIZE);

		/* plane 1 */
		if (itr >= cx1-len1/2 && itr <= cx1+len1/2) {
		    ++tfe1;

		    /* fit is fractional sample of plane intersection */
		    fit = ct1 - ( cx1 - itr ) * dip1 / msps; 
		    if (fit >= 0 && fit <= (float) nt) {

			/* linear interpolation */
			itless = fit;
			eps = fit - itless;
			itmore = fit + 1;
			tr.data[itless] += 1.0 - eps;	 
			tr.data[itmore] += eps;	 

			/* taper option */
			if (taper == 1) {
			  /* last point */
			  if (tfe1 == 1 || tfe1 == len1 + 1) {
				tr.data[itless] /= 6.0;	 
				tr.data[itmore] /= 6.0;	 
			  } 
			  /*  next-to-last point */
			  if (tfe1 == 2 || tfe1 == len1) {
				tr.data[itless] /= 3.0;	 
				tr.data[itmore] /= 3.0;	 
			  } 
		    }
		  }
		}

		/*  plane 2  */
		if (npl > 1) {
		  if (itr >= cx2-len2/2 && itr <= cx2+len2/2) {
		    ++tfe2;

		    /* fit is fractional sample of plane intersection */
		    fit = ct2 - ( cx2 - itr ) * dip2 / msps; 
		    if (fit >= 0 && fit <= (float) nt) {

			/* linear interpolation */
			itless = fit;
			eps = fit - itless;
			itmore = fit + 1;
			tr.data[itless] += 1.0 - eps;	 
			tr.data[itmore] += eps;	 

			/* taper option */
			if (taper == 1) {
			  /* last point */
			  if (tfe2 == 1 || tfe2 == len2 + 1) {
				tr.data[itless] /= 6.0;	 
				tr.data[itmore] /= 6.0;	 
			  } 
			  /*  next-to-last point */
			  if (tfe2 == 2 || tfe2 == len2) {
				tr.data[itless] /= 3.0;	 
				tr.data[itmore] /= 3.0;	 
			  } 
		        }
		    }
		  }
		}

		/* plane 3  */
		if (npl > 2) {
		  if (itr >= cx3-len3/2 && itr <= cx3+len3/2) {
		    ++tfe3;

		    /* fit is fractional sample of plane intersection */
		    fit = ct3 - ( cx3 - itr ) * dip3 / msps; 
		    if (fit >= 0 && fit <= (float) nt) {

			/* linear interpolation */
			itless = fit;
			eps = fit - itless;
			itmore = fit + 1;
			tr.data[itless] += 1.0 - eps;	 
			tr.data[itmore] += eps;	 

			/* taper option */
			if (taper == 1) {
			  /* last point */
			  if (tfe3 == 1 || tfe3 == len3 + 1) {
				tr.data[itless] /= 6.0;	 
				tr.data[itmore] /= 6.0;	 
			  } 
			  /* next-to-last point */
			  if (tfe3 == 2 || tfe3 == len3) {
				tr.data[itless] /= 3.0;	 
				tr.data[itmore] /= 3.0;	 
			  } 
		        }
		    }
		  }
		}

		/* fill tracl/tracr headers and put trace out */
		tr.tracl = tr.tracr = itr + 1;
		su2cs->putTrace(&tr);
	}
	
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
