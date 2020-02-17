/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUINTVEL: $Revision: 1.14 $ ; $Date: 2011/11/16 17:43:20 $		*/

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
std::string sdoc_suintvel =
"									"
" SUINTVEL - convert stacking velocity model to interval velocity model	"
"									"
" suintvel vs= t0= outpar=/dev/tty					"
"									"
" Required parameters:					        	"
"	vs=	stacking velocities 					"
"	t0=	normal incidence times		 			"
"									"
" Optional parameters:							"
"	mode=0			output h= v= ; =1 output v=  t= 	"
"	outpar=/dev/tty		output parameter file in the form:	"
"				h=layer thicknesses vector		"
"				v=interval velocities vector		"
"				....or ...				"
"				t=vector of times from t0		"
"				v=interval velocities vector		"
"									"
" Examples:								"
"    suintvel vs=5000,5523,6339,7264 t0=.4,.8,1.125,1.425 outpar=intpar	"
"									"
"    suintvel par=stkpar outpar=intpar					"
"									"
" If the file, stkpar, contains:					"
"    vs=5000,5523,6339,7264						"
"    t0=.4,.8,1.125,1.425						"
" then the two examples are equivalent.					"
"									"
" Note: suintvel does not have standard su syntax since it does not	"
"      operate on seismic data.  Hence stdin and stdout are not used.	"
"									"
" Note: may go away in favor of par program, velconv, by Dave		"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suintvel {


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
 *	Define:
 *		t0by2[k] = 0.5 * t0[k] = Sum h[i]/v[i]
 *		vh[k] = vs[k]*vs[k]*t0by2[k] = Sum v[i]*h[i]
 *	Then:
 *		dt[i] = h[i]/v[i] = t0by2[i] - t0by2[i-1]
 *		dvh[i] = h[i]*v[i] = vh[i] - vh[i-1]
 *		h[i] = sqrt(dvh[i] * dt[i])
 *		v[i] = sqrt(dvh[i] / dt[i])
 *
 *
 */
/**************** end self doc *******************************************/

void* main_suintvel( void* args )
{
	register float *v=NULL;		/* interval velocities		*/
	register float *h=NULL;		/* layer thicknesses at the cmp	*/
	register float *vs=NULL;	/* stacking velocities		*/
	register float *t0=NULL;	/* zero incidence times		*/
	register int i;		/* counter				*/
	int n;			/* number of layers			*/
	float t1, t2;		/* temporaries for one-way times	*/
	float v1, v2;		/* temporaries for stacking v's		*/
	float dt;		/* temporary for t0/2 difference	*/
	float dvh;		/* temporary for v*h difference		*/
	cwp_String outpar;	/* name of file holding output parfile	*/
	FILE *outparfp;		/* ... its file pointer			*/
	int mode;		/* mode=0  h= v= ; mode=1 t= v= 	*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suintvel );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	outpar = "/dev/tty" ;	parObj.getparstring("outpar", &outpar);
	outparfp = efopen(outpar, "w");


	/* Allocate space for the model */
	if ((n = parObj.countparval("t0"))) {
		v  = ealloc1float(n);
		h  = ealloc1float(n);
		vs = ealloc1float(n);
		t0 = ealloc1float(n);
	} else throw cseis_geolib::csException("no t0's specified");

	/* Get the normal incidence times and stacking velocities */
	if (n != parObj.getparfloat("t0", t0))
		throw cseis_geolib::csException("expected %d intervals", n);
	if (n != parObj.getparfloat("vs", vs))
		throw cseis_geolib::csException("expected %d velocities", n);
	if (!parObj.getparint("mode", &mode))		mode = 0;
        parObj.checkpars();

	/* Check that vs's and t0's are positive */
	for (i = 0; i < n; i++) {
		if (vs[i] <= 0.0)
			throw cseis_geolib::csException("vs's must be positive: vs[%d] = %f", i, vs[i]);
		if (t0[i] < 0.0)
			throw cseis_geolib::csException("t0's must be positive: t0[%d] = %f", i, t0[i]);
	}

	/* Compute h(i), v(i) */
	h[0] = 0.5 * vs[0] * t0[0];
	v[0] = vs[0];
	for (i = 1; i < n; i++) {
		t2 = 0.5 * t0[i]; t1 = 0.5 * t0[i-1];
		v2 = vs[i]; v1 = vs[i-1];
		dt = t2 - t1;
		dvh = v2*v2*t2 - v1*v1*t1;
		h[i] = sqrt(dvh * dt);
		v[i] = sqrt(dvh / dt);
	}

	/* Make par file */
	if (!mode) {
		fprintf(outparfp, "h=");
		for (i = 0; i < n - 1; i++) {
			fprintf(outparfp, "%g,", h[i]);
		}
		fprintf(outparfp, "%g\n", h[n-1]);
	} else {
		fprintf(outparfp, "t=");
		for (i = 0; i < n - 1; i++) {
			fprintf(outparfp, "%g,", t0[i]);
		}
		fprintf(outparfp, "%g\n", t0[n-1]);
	}
	fprintf(outparfp, "v=");
	for (i = 0; i < n - 1; i++) {
		fprintf(outparfp, "%g,", v[i]);
	}
	fprintf(outparfp, "%g\n", v[n-1]);


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
