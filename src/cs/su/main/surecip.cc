/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SURECIP: $Revision: 1.14 $ ; $Date: 2011/11/16 23:14:54 $	*/

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
std::string sdoc_surecip =
"                                				"
" SURECIP - sum opposing offsets in prepared data (see below)	"
"                                				"
" surecip <stdin >stdout	 		               	"
"                                				"
" Sum traces with equal positive and negative offsets (i.e. assume"
" reciprocity holds). 						"
"                                				"
" Usage:							"
"	suabshw <data >absdata					"
"	susort cdp offset <absdata | surecip >sumdata		"
"                                				"
" Note that this processing stream can be simply evoked by:	"
"                                				"
"	recip data sumdata					"
"                                				"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace surecip {


/* Credits:
 *	SEP: Shuki Ronen
 *	CWP: Jack Cohen
 *
 * Caveat:
 *	The assumption is that this operation is not a mainstay processing
 *	item.  Hence the recommended implemention via the 'recip' shell
 *	script.  If it becomes a mainstay, then a much faster code can
 *	quickly drummed up by incorporating portions of suabshw and
 *	susort.
 *
 * Trace header fields accessed: ns
 * Trace header fields modified: nhs, tracl, sx, gx
 */
/**************** end self doc ***********************************/


segy intrace, outtrace;

void* main_surecip( void* args )
{
	int cdpindex;	/* index of cdp header word		*/
	int offindex;	/* index of offset header word		*/
	int nt;		/* number of data points on trace	*/
	int nsegy;	/* number of bytes in the segy		*/
	Value cdpval;	/* value of cdp in current gather	*/
	Value offval;	/*  ... same for offset			*/
	Value cdpvalnew;/* value of cdp in trace being treated	*/
	Value offvalnew;/* ... same for offset			*/
	int newtracl;	/* tracl for stacked traces		*/
	int fold;	/* number of traces with same offset 	*/
	float ffold;	/* ... cast to float			*/
	int norm;	/* norm=1 => divide by fold		*/
	int itmp;	/* temporary for swap of sx, gx keys	*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_surecip );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Set parameter */
	if (!parObj.getparint("norm", &norm))	norm = 1;

        parObj.checkpars();

	/* Get indices */
	cdpindex = getindex("cdp");
	offindex = getindex("offset");

	/* Set up for first trace */
	nsegy = cs2su->getTrace(&intrace);
	nt = intrace.ns;
	gethval(&intrace, cdpindex, &cdpval);
	gethval(&intrace, offindex, &offval);
	memcpy( (void *) &outtrace, (const void *) &intrace, nsegy);

	newtracl = 1;		/* Global initialization */
	fold = 1;		/* Will be re-initialized for each gather */

	/* Loop over traces */
	while (nsegy) {		     /* While previous trace non-empty */
		nsegy = cs2su->getTrace(&intrace);
		gethval(&intrace, cdpindex, &cdpvalnew);
		gethval(&intrace, offindex, &offvalnew);
		if (valcmp("l", cdpval, cdpvalnew) || !nsegy) {	
			/* Either cdpval and cdpvalnew differ,       */
			/* indicating a new gather or nsegy is zero, */
			/* indicating the end of the traces.         */

			/* Add header info and output leftover stack */
			outtrace.nhs = fold;
			outtrace.tracl = newtracl++;
			if (outtrace.sx > outtrace.gx) {
				itmp = outtrace.sx;
				outtrace.sx = outtrace.gx;
				outtrace.gx = itmp;
			}
			if (norm) {
				ffold = (float) fold;
				if (fold != 1) {
				    register int i;
				    for (i = 0; i < nt; ++i)
					outtrace.data[i] /= ffold;
				}
			}
			su2cs->putTrace(&outtrace);
	
			/* Set up for next gather */
			memcpy( (void *) &outtrace,
					(const void *) &intrace, nsegy);
			fold = 1;
			cdpval = cdpvalnew;
			offval = offvalnew;

		} else {	/* still in same cdp gather */
			if (valcmp("l", offval, offvalnew)) {
				/* offval and offvalnew differ */

				/* Add header info and output stack */
				outtrace.nhs = fold;
				outtrace.tracl = newtracl++;
				if (outtrace.sx > outtrace.gx) {
					itmp = outtrace.sx;
					outtrace.sx = outtrace.gx;
					outtrace.gx = itmp;
				}
				if (norm) {
					ffold = (float) fold;
					if (fold != 1) {
					    register int i;
					    for (i = 0; i < nt; ++i)
						outtrace.data[i] /= ffold;
					}
				}
				su2cs->putTrace(&outtrace);

				/* Set up for next offset */
				memcpy( (void *) &outtrace,
				       (const void *) &intrace, nsegy);
				fold = 1;
				offval = offvalnew;

			} else { /* same offset within this cdp */

				register int i;
				for (i = 0; i < nt; ++i)
					outtrace.data[i] += intrace.data[i];

				fold++;
			}
		}
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
