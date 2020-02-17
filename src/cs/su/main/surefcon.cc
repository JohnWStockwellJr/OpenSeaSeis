/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUREFCON: $Revision: 1.6 $ ; $Date: 2011/11/16 17:37:27 $	*/

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
/*********************** self documentation **************************/
std::string sdoc_surefcon =
"									"
" SUREFCON -  Convolution of user-supplied Forward and Reverse		"
"		refraction shots using XY trace offset in reverse shot	"
"									"
"	surefcon <forshot sufile=revshot  xy=trace offseted  >stdout	"
"									"
" Required parameters:						 	"
" sufile=	file containing SU trace to use as reverse shot		"
" xy=		Number of traces offseted from the 1st trace in sufile	"
"									"
" Optional parameters:						 	"
" none								 	"
"									"
" Trace header fields accessed: ns					"
" Trace header fields modified: ns					"
"									"
" Notes:								"
" This code implements the Refraction Convolution Section (RCS)	method	"
" of generalized reciprocal refraction traveltime analysis developed by "
" Derecke Palmer and Leoni Jones.					"
"									"
" The time sampling interval on the output traces is half of that on the"
" traces in the input files.		  	"
"									"
" Example:								"
"									"
"	 surefcon <DATA sufile=DATA xy=1 | ...				"
"									"
" Here, the su data file, \"DATA\" convolved the nth trace by		"
" (n+xy)th trace in the same file					"
"									"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace surefcon {


/* Credits: (based on suconv)
 *	CWP: Jack K. Cohen, Michel Dietrich
 *	UNSW: D. Palmer, K.T. LEE
 *  CAVEATS: no space-variable or time-variable capacity.
 *	The more than one trace allowed in sufile is the
 *	beginning of a hook to handle the spatially variant case.
 *
 * Trace header fields accessed: ns
 * Trace header fields modified: ns
 * Notes:
 * This code implements the refraction convolution 
 * section (RCS) method
 * method described in:
 *
 * Palmer, D, 2001a, Imaging refractors with the convolution section,
 *           Geophysics 66, 1582-1589.
 * Palmer, D, 2001b, Resolving refractor ambiguities with amplitudes,
 *           Geophysics 66, 1590-1593.
 *
 * Exploration Geophysics (2005) 36, 18­25
 * Butsuri-Tansa (Vol. 58, No.1)
 * Mulli-Tamsa (Vol. 8,
 *    A simple approach to refraction statics with the 
 * Generalized Main Reciprocal Method and the Refraction 
 * Convolution Section Heading
 *        by Derecke Palmer  Leonie Jones
 *
 */
/**************** end self doc ********************************/

segy intrace, outtrace, sutrace;

void* main_surefcon( void* args )
{
	int nt;		/* number of points on input traces	*/
	int ntout;	/* number of points on output traces	*/
	int xy;		/* the offset number for GRM		*/
	float *forshot;	/* forward shot			*/
	int nforshot;	/* length of input wavelet in samples  */
	cwp_String sufile;	/* name of file of forward SU traces	*/
	FILE *fp;		/* ... its file pointer		*/
	int itr;		/* trace counter			 */

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_surefcon );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get info from first trace */ 
	if (!cs2su->getTrace(&intrace) ) throw cseis_geolib::csException("can't get 1st reverse shot trace");
	nt = intrace.ns;

	/* Default parameters;  User-defined overrides */
	if (!parObj.getparint("xy", &xy) )  xy = 0;
		
	/* Check xy values */
	if (xy < 0)	throw cseis_geolib::csException("xy=%d should be positive", xy);
	 
	if (!parObj.getparstring("sufile", &sufile)) {
		throw cseis_geolib::csException("must specify sufile= desired forward shot");
	} else {
                parObj.checkpars();

		/* Get parameters and set up forshot array */
		fp = efopen(sufile, "r");
	 	for (itr = 0; itr <= xy; ++itr) {
			if (!fgettr(fp, &sutrace) ) { 
				throw cseis_geolib::csException("can't get 1st requested forward trace");
			}
	 	}	

	 	nforshot = sutrace.ns;
	 	forshot = ealloc1float(nforshot);

		 /* Set output trace length */
		 ntout = nt + nforshot - 1;
		
		 /* Main loop over reverse shot traces */
		 do {  
			warn("rev==%d\t , for=%d", intrace.tracf, sutrace.tracf);
			 memcpy((void *) forshot,
			(const void *) sutrace.data, nforshot*FSIZE);
	
			/* Convolve for shot with reverse shot trace */
			convolve_cwp(nforshot, 0, forshot,
					 nt, 0, intrace.data, 
					ntout, 0, outtrace.data);	

			/* Output convolved trace */
			memcpy((void *) &outtrace, (const void *) &intrace, HDRBYTES);
			outtrace.ns = ntout;
			outtrace.dt = outtrace.dt/2;
			/*outtrace.cdp = 2*intrace.tracf + xy;*/ 
			fprintf(stderr,"out_cdp=%d\n", 2*intrace.tracf + xy);
			su2cs->putTrace(&outtrace);
		
		 } while ( cs2su->getTrace(&intrace) && fgettr(fp, &sutrace) ); 
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
