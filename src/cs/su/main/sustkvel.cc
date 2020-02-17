/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSTKVEL: $Revision: 1.14 $ ; $Date: 2011/11/16 17:43:20 $		*/

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

/*********************** self documentation ******************************/
std::string sdoc_sustkvel =
"									"
" SUSTKVEL - convert constant dip layer interval velocity model to the	"
"	   stacking velocity model required by sunmo			"
"									"
" sustkvel v= h= dip=0.0 outpar=/dev/tty				"
"									"
" Required parameters:					        	"
"	v=	interval velocities 					"
"	h=	layer thicknesses at the cmp	 			"
"									"
" Optional parameters:							"
"	dip=0.0			(constant) dip of the layers (degrees)	"
"	outpar=/dev/tty		output parameter file in the form	"
"				required by sunmo:			"
"				tv=zero incidence time pick vector	"
"				v=stacking velocities vector		"
"									"
" Examples:								"
"    sustkvel v=5000,6000,8000,10000 h=1000,1200,1300,1500 outpar=stkpar"
"    sunmo <data.cdp par=stkpar >data.nmo				"
"									"
"    sustkvel par=intpar outpar=stkpar					"
"    sunmo <data.cdp par=stkpar >data.nmo				"
"									"
" If the file, intpar, contains:					"
"    v=5000,6000,8000,10000						"
"    h=1000,1200,1300,1500						"
" then the two examples are equivalent.  The created parameter file,	"
" stkpar, is in the form of the velocity model required by sunmo.	"
"									"
" Note: sustkvel does not have standard su syntax since it does not	"
"      operate on seismic data.  Hence stdin and stdout are not used.	"
"									"
" Caveat: Does not accept a series of interval velocity models to	"
"	produce a variable velocity file for sunmo.			"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sustkvel {


/* Credits:
 *	CWP: Jack 
 *
 * Technical Reference:
 *	The Common Depth Point Stack
 *	William A. Schneider
 *	Proc. IEEE, v. 72, n. 10, p. 1238-1254
 *	1984
 *
 * Formulas:
 *    	Note: All sums on i are from 1 to k
 *
 *	From Schneider:
 *	Let h[i] be the ith layer thickness measured at the cmp and
 *	v[i] the ith interval velocity.
 *	Set:
 *		t[i] = h[i]/v[i]
 *		t0[k] = 2 Sum t[i] * cos(dip)
 *		vs[k] = (1.0/cos(dip)) sqrt(Sum v[i]*v[i]*t[i] / Sum t[i])
 *	Define:
 *		t0by2[k] = Sum h[i]/v[i]
 *		vh[k]    = Sum v[i]*h[i]
 *	Then:
 *		t0[k] = 2 * t0by2[k] * cos(dip)
 *		vs[k] = sqrt(vh[k] / t0by2[k]) / cos(dip)
 *
 *
 */
/**************** end self doc *******************************************/

void* main_sustkvel( void* args )
{
	register float *v=NULL;		/* interval velocities		*/
	register float *h=NULL;		/* layer thicknesses at the cmp	*/
	register float *vs=NULL;	/* stacking velocities		*/
	register float *t0=NULL;	/* zero incidence times		*/
	register int i;		/* counter				*/
	int n;			/* number of layers			*/
	float t0by2;		/* temporary for t0/2 sum		*/
	float vh;		/* temporary for vh sum			*/
	float dip;		/* (constant) layer dip (degrees)	*/
	float cosdip;		/* cos(dip)				*/
	float ocosdip;		/* 1.0/cos(dip)				*/
	cwp_String outpar;	/* name of file holding output parfile	*/
	FILE *outparfp;		/* ... its file pointer			*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sustkvel );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get parameters */
	if (!parObj.getparfloat("dip", &dip))	dip = 0.0;
	dip = dip * PI / 180.0;
	cosdip = cos(dip);
	ocosdip = 1.0 / cosdip;

	if (!parObj.getparstring("outpar", &outpar))	outpar = "/dev/tty";
	outparfp = efopen(outpar, "w");


	/* Allocate space for the model */
	if ((n = parObj.countparval("v"))) {
		v  = ealloc1float(n);
		h  = ealloc1float(n);
		vs = ealloc1float(n);
		t0 = ealloc1float(n);
	} else throw cseis_geolib::csException("no v's specified");

	/* Get the intervals and interval velocities */
	if (n != parObj.getparfloat("h", h))
		throw cseis_geolib::csException("expected %d intervals", n);
	if (n != parObj.getparfloat("v", v))
		throw cseis_geolib::csException("expected %d velocities", n);

        parObj.checkpars();

	/* Check that v's and h's are positive */
	for (i = 0; i < n; i++) {
		if (v[i] <= 0.0)
			throw cseis_geolib::csException("v's must be positive: v[%d] = %f", i, v[i]);
		if (h[i] <= 0.0)
			throw cseis_geolib::csException("h's must be positive: h[%d] = %f", i, h[i]);
	}

	/* Compute T0(i), Vs(i) */
	t0by2 = 0.0;
	vh = 0.0;
	for (i = 0; i < n; i++) {
		t0by2 += h[i]/v[i];
		vh += v[i] * h[i];
		t0[i] = 2.0 * t0by2 * cosdip;
		vs[i] = sqrt(vh/t0by2) * ocosdip;
	}

	/* Make par file */
	fprintf(outparfp, "tnmo=");
	for (i = 0; i < n - 1; i++) {
		fprintf(outparfp, "%g,", t0[i]);
	}
	fprintf(outparfp, "%g\n", t0[n-1]);

	fprintf(outparfp, "vnmo=");
	for (i = 0; i < n - 1; i++) {
		fprintf(outparfp, "%g,", vs[i]);
	}
	fprintf(outparfp, "%g\n", vs[n-1]);


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
