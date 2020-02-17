/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUCONV: $Revision: 1.18 $ ; $Date: 2011/11/16 17:37:27 $		*/

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
std::string sdoc_suconv =
" 									"
" SUCONV - convolution with user-supplied filter			"
" 									"
" suconv <stdin >stdout  filter= [optional parameters]			"
" 									"
" Required parameters: ONE of						"
" sufile=		file containing SU trace to use as filter	"
" filter=		user-supplied convolution filter (ascii)	"
" 									"
" Optional parameters:							"
" panel=0		use only the first trace of sufile		"
" 			=1 convolve corresponding trace in sufile with	"
" 			trace in input data				"
" 									"
" Trace header fields accessed: ns					"
" Trace header fields modified: ns					"
" 									"
" Notes: It is quietly assumed that the time sampling interval on the	"
" single trace and the output traces is the same as that on the traces	"
" in the input file.  The sufile may actually have more than one trace,	"
" but only the first trace is used in panel=0. In panel=1 the corresponding"
" trace from the sufile are convolved with its counterpart in the data.	"
" Caveat, in panel=1 there have to be at least as many traces in sufile	"
" as in the input data. If not, a warning is returned, and later traces	"
" in the dataset are returned unchanged.				"
" 									"
" Examples:								"
"	suplane | suwind min=12 max=12 >TRACE				"
"	suconv<DATA sufile=TRACE | ...					"
" Here, the su data file, \"DATA\" is convolved trace by trace with the"
" the single su trace, \"TRACE\".					"
" 									"
"	suconv<DATA filter=1,2,1 | ...					"
" Here, the su data file, \"DATA\" is convolved trace by trace with the"
" the filter shown.							"
" 									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suconv {


/* Credits:
 *	CWP: Jack K. Cohen, Michel Dietrich
 *
 *  CAVEATS: no space-variable or time-variable capacity.
 *     The more than one trace allowed in sufile is the
 *     beginning of a hook to handle the spatially variant case.
 *
 * Trace header fields accessed: ns
 * Trace header fields modified: ns
 */
/**************** end self doc *******************************************/

segy intrace, outtrace, sutrace;

void* main_suconv( void* args )
{
	int nt;			/* number of points on input traces	*/
	int ntout;		/* number of points on output traces	*/
	float *filter=NULL;	/* filter coefficients			*/
	int nfilter=0;		/* length of input wavelet in samples	*/
	cwp_String sufile="";	/* name of file containing one SU trace */
	FILE *fp=NULL;		/* ... its file pointer			*/
	int delrtf=0;		/* delrt from traces in sufile		*/

	int panel;		/* operational panel of the program	*/
	cwp_Bool issufile=cwp_false;	/* is sufile set?		*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suconv );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	if (!parObj.getparint("panel",&panel))		panel=0;

	/* Get info from first trace */ 
	if (!cs2su->getTrace(&intrace) ) throw cseis_geolib::csException("can't get first trace");
	nt = intrace.ns;

	/* Get parameters and set up filter array */
	if (!parObj.getparstring("sufile", &sufile) ) {
		if (!(nfilter = parObj.countparval("filter")))
			throw cseis_geolib::csException("must specify filter= desired filter");
		filter = ealloc1float(nfilter);	parObj.getparfloat("filter", filter);
	} else if (panel==0) { /* if panel=0 only use the first trace */
		fp = efopen(sufile, "r");
		fgettr(fp, &sutrace);

                /* HD: get delrt from filter sufile */
                delrtf = sutrace.delrt;

		nfilter = sutrace.ns;
		filter = ealloc1float(nfilter);
		memcpy((void *) filter,
			(const void *) sutrace.data, nfilter*FSIZE);
		issufile=cwp_true; 
	} else { /* if panel=1 use each trace */
		fp = efopen(sufile, "r");
		issufile=cwp_true; 
	}
        parObj.checkpars();

	/* Set output trace length */
	ntout = nt + nfilter - 1;

	/* Main loop over traces */
	do {
		/* if panel=1 and sufile is defined */
		if ((panel==1) && (issufile==cwp_true)) {
			fgettr(fp, &sutrace);

                	/* HD: get delrt from filter sufile */
                	delrtf = sutrace.delrt;

			nfilter = sutrace.ns;
			filter = ealloc1float(nfilter);
			memcpy((void *) filter,
				(const void *) sutrace.data, nfilter*FSIZE);
			ntout = nt + nfilter - 1;
		}

		/* Convolve filter with trace */
		convolve_cwp(nfilter, 0, filter,
		     nt, 0, intrace.data, 
                     ntout, 0, outtrace.data);        


		/* Output filtered trace */
		memcpy((void *) &outtrace, (const void *) &intrace, HDRBYTES);
		outtrace.ns = ntout; 

		/* HD: update delrt */
		outtrace.delrt = intrace.delrt+delrtf;

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
