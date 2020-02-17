/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSTRIP: $Revision: 1.20 $ ; $Date: 2011/11/16 22:10:29 $	*/

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

/*********************** self documentation **********************/
std::string sdoc_sustrip =
" 								"
" SUSTRIP - remove the SEGY headers from the traces		"
" 								"
" sustrip <stdin >stdout head=/dev/null outpar=/dev/tty ftn=0	"
" 								"
" Required parameters:						"
" 	none							"
" 								"
" Optional parameters:						"
" 	head=/dev/null		file to save headers in		"
" 								"
" 	outpar=/dev/tty		output parameter file, contains:"
" 				number of samples (n1=)		"
" 				number of traces (n2=)		"
" 				sample rate in seconds (d1=)	"
" 								"
" 	ftn=0			Fortran flag			"
" 				0 = write unformatted for C	"
" 				1 = ... for Fortran		"
" 								"
" Notes:							"
" Invoking head=filename will write the trace headers into filename."
" You may paste the headers back onto the traces with supaste	"
" See:  sudoc  supaste 	 for more information 			"
" Related programs: supaste, suaddhead				"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sustrip {


/* Credits:
 *	SEP: Einar Kjartansson  c. 1985
 *	CWP: Jack K. Cohen        April 1990
 *
 * Trace header fields accessed: ns, dt
 */
/**************** end self doc ***********************************/


segy tr;

void* main_sustrip( void* args )
{
	cwp_String head;	/* name of file holding headers		*/
	FILE *headfp; 		/* ... its file pointer			*/
	cwp_String outpar;	/* name of file holding output parfile	*/
	FILE *outparfp;		/* ... its file pointer			*/
	int ns;			/* number of data samples on the segys	*/
	size_t nsbytes;		/* ... in bytes				*/
	int ftn;		/* fortran flag				*/
	int ntr = 0;		/* number of traces written		*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sustrip );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	switch(filestat(STDOUT)) {
	case BADFILETYPE:
		warn("stdout is illegal filetype");
		pagedoc();
	break;
	case TTY:
		warn("stdout can't be tty");
		pagedoc();
	break;
	default: /* other cases are OK */
	break;
	}

	/* Get parameters */
	if (!parObj.getparstring("head"  , &head))	head   = "/dev/null";
	if (!parObj.getparstring("outpar", &outpar))	outpar = "/dev/tty" ;
	if (!parObj.getparint   ("ftn"   , &ftn))	ftn = 0;
	if (ftn != 0 && ftn != 1)  throw cseis_geolib::csException("ftn=%d must be 0 or 1", ftn);


        parObj.checkpars();
	/* Open files to save headers and parameters */
	headfp = efopen(head, "w");
	outparfp = efopen(outpar, "w");


	/* Get info from first trace */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	ns = tr.ns;
	nsbytes = ns * FSIZE;


	/* Write the data portion of the records--if the ftn	*/
	/* option is selected, write an int before and after	*/
	/* each trace giving the length of the trace in bytes	*/
	/* as per the Fortran unformatted record format.	*/
	do {

		if (ftn) efwrite(&nsbytes, ISIZE, 1, stdout);
		switch(tr.trid) {
		case CHARPACK:  efwrite(tr.data, sizeof(char), ns, stdout);
		break;
		case SHORTPACK: efwrite(tr.data, sizeof(short), ns, stdout);
		break;
		default:        efwrite(tr.data, FSIZE, ns, stdout);
		}     
		if (ftn) efwrite(&nsbytes, ISIZE, 1, stdout);

		efwrite(&tr, 1, HDRBYTES, headfp);

		++ntr;

	} while (cs2su->getTrace(&tr));

	/* Make par file for headerless file */
	fprintf(outparfp, "n1=%d n2=%d d1=%f\nnt=%d ntr=%d dt=%f\nns=%d\n",
		tr.ns, ntr, ((double) tr.dt)/1000000.0,
		tr.ns, ntr, ((double) tr.dt)/1000000.0,
		tr.ns);


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
