/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUPACK2: $Revision: 1.22 $ ; $Date: 2011/11/16 17:38:58 $	*/

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
std::string sdoc_supack2 =
"								"
" SUPACK2 - pack segy trace data into 2 byte shorts		"
"								"
" supack2 <segy_file >packed_file	gpow=0.5 		"
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
namespace supack2 {


/* Credits:
 *	CWP: Jack K. Cohen, Shuki Ronen, Brian Sumner
 *
 * Revised: 7/4/95  Stewart A. Levin  Mobil
 *          Changed encoding to ensure 2 byte length (short is
 *	    8 bytes on Cray).
 *
 * Caveats:
 *	This program is for single site use.  Use segywrite to make
 *	a portable tape.
 *
 *	We are storing the local header words, ungpow and unscale,
 *	required by suunpack2 as floats.
 *	
 * Notes:
 *	ungpow and unscale are defined in segy.h
 *	trid = SHORTPACK is defined in su.h and segy.h
 *
 * Trace header fields accessed: ns
 * Trace header fields modified: ungpow, unscale, trid
 */
/**************** end self doc ***********************************/


#define GPOW	0.5	/* default power parameter */

segy tr;	/* on  input: SEGY hdr & (float) trace data */
		/* on output: data as 2-byte shorts          */

void* main_supack2( void* args )
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
	su2cs->setSUDoc( sdoc_supack2 );
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
		 * Since the shorts take less room than the floats,
		 * we don't overwrite.
		 *
		 * Note that the segy field tr.data is declared as
		 * floats, so we need to invent a pointer for the
		 * short array which is actually there. */

		register unsigned char *otr = (unsigned char *) tr.data;
		register int i,j;
		register signed   int si;
		register unsigned int ui;
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
					    pow(val, gpow) : -pow(-val, gpow);
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
		scale = absmax ? SHRT_MAX/absmax : 0.0;
		tr.unscale = absmax ? 1.0/scale : 0.0;

		/* Apply the scale and load in short data
		 * Note: the essence of the code below is:
		 * for (i = 0; i < nt; ++i) { 
		 *	tr.data[i] *= scale;
		 *      otr[i] = (short) tr.data[i];
		 * }
		 * but this assumes shorts are 2 bytes, so isn't portable */
		for (i = 0, j=0; i < nt; ++i) { 
			tr.data[i] *= scale;
			si = (signed int) tr.data[i];
			ui = (si>>8)&255;
			otr[j++] = (unsigned char) ui;
			ui = si&255;
			otr[j++] = (unsigned char) ui;
		}

		/* Write trace ID as the packed short code number */
		tr.trid = SHORTPACK;

		/* Output the "segy" with shorts in the data array */
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
