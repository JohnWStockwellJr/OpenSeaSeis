/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUPACK1: $Revision: 1.20 $ ; $Date: 2011/11/16 17:38:58 $	*/

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
std::string sdoc_supack1 =
"								"
" SUPACK1 - pack segy trace data into chars			"
"								"
" supack1 <segy_file >packed_file	gpow=0.5 		"
"								"
" Required parameters:						"
"	none							"
"						        	"
" Optional parameter: 						"
"	gpow=0.5	exponent used to compress the dynamic	"
"			range of the traces			"
"								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace supack1 {


/* Credits:
 *	CWP: Jack K. Cohen, Shuki Ronen, Brian Sumner
 *
 * Caveats:
 *	This program is for single site use.  Use segywrite to make
 *	a portable tape.
 *
 *	We are storing the local header words, ungpow and unscale,
 *	required by suunpack1 as floats.  Although not essential
 *	(compare the handling of such fields as dt), it allows us
 *	to demonstrate the convenience of using the natural data type.
 *	In any case, the data itself is non-portable floats in general,
 *	so we aren't giving up any intrinsic portability.
 *	
 * Notes:
 *	ungpow and unscale are defined in segy.h
 *	trid = CHARPACK is defined in su.h and segy.h
 *
 * Trace header fields accessed: ns
 * Trace header fields modified: ungpow, unscale, trid
 */
/**************** end self doc ***********************************/


#define GPOW	0.5	/* default power parameter */

segy tr;	/* on  input: SEGY hdr & (float) trace data */
		/* on output: data as signed chars          */

void* main_supack1( void* args )
{
	float gpow;
	int nt;
	cwp_Bool isone, ishalf;
	float f_one = 1.0;
	float f_half = 0.5;


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_supack1 );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get parameters */
	if (!parObj.getparfloat("gpow", &gpow)) gpow = GPOW;
        parObj.checkpars();
	if (gpow <= 0.0) throw cseis_geolib::csException("gpow = %g must be positive", gpow);
	isone  = CLOSETO(gpow, f_one);
	ishalf = CLOSETO(gpow, f_half);

	/* Get number of time samples from first trace */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;

	/* Main loop over segy traces */
	do {
		/* Point output trace at the trace data and pack.
		   Since the chars take less room than the floats,
		   we don't overwrite.
	
		   Note that the segy field tr.data is declared as
		   floats, so we need to invent a pointer for the
		   char array which is actually there.
		*/

		register signed char *otr = (signed char *) tr.data;
		register int i;
		register float absmax;
		register float scale;

		/* Power transform to decrease dynamic range */
		if (!isone) {
			register float val;

			if (ishalf) {
				for (i = 0; i < nt; ++i) {
					val = tr.data[i];
					tr.data[i] = (val >= 0.0) ?
						sqrt(val) : -sqrt(-val);
				}
			} else {
				for (i = 0; i < nt; ++i) {
					val = tr.data[i];
					tr.data[i] = (val >= 0.0) ?
						pow(val, gpow) :
							-pow(-val, gpow);
				}
			}
		}

		/* Store "ungpow" factor */
		tr.ungpow = 1.0/gpow;

		/* Read trace data and get absmax */
		absmax = ABS(tr.data[0]);
		for (i = 1; i < nt; ++i)
			absmax = MAX(absmax, ABS(tr.data[i]));

		/* Compute scale factor and store "unscale" factor */
		/* If max is zero, then put scale and unscale to zero too */
		scale = absmax ? SCHAR_MAX/absmax : 0.0;
		tr.unscale = absmax ? 1.0/scale : 0.0;

		/* Apply the scale and load in char data */
		for (i = 0; i < nt; ++i) { 
			tr.data[i] *= scale;
			otr[i] = (signed char) tr.data[i];
		}

		/* Write trace ID as the packed char code number */
		tr.trid = CHARPACK;

		/* Output the "segy" with chars in the data array */
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
