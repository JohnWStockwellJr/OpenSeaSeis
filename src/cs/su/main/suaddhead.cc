/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUADDHEAD: $Revision: 1.21 $ ; $Date: 2011/11/16 22:10:29 $		*/

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
std::string sdoc_suaddhead =
" 									"
" SUADDHEAD - put headers on bare traces and set the tracl and ns fields"
" 									"
" suaddhead <stdin >stdout ns= ftn=0					"
" 									"
" Required parameter:							"
" 	ns=the number of samples per trace				"
" 									"
" Optional parameter:							"
#ifdef SU_LINE_HEADER
"	head=           file to read headers in				"
"                       not supplied --  will generate headers 		"
"                       given        --  will read in headers and attach"
"                                        floating point arrays to form 	"
"                                        traces 			" 
"                       (head can be created via sustrip program)	"
#endif
" 	ftn=0		Fortran flag					"
" 			0 = data written unformatted from C		"
" 			1 = data written unformatted from Fortran	"
"       tsort=3         trace sorting code:				"
"                                1 = as recorded (no sorting)		"
"                                2 = CDP ensemble			"
"                                3 = single fold continuous profile	"
"                                4 = horizontally stacked		" 
"       ntrpr=1         number of data traces per record		"
"                       if tsort=2, this is the number of traces per cdp" 
" 									"
" Trace header fields set: ns, tracl					"
" Use sushw/suchw to set other needed fields.				"
" 									"
" Caution: An incorrect ns field will munge subsequent processing.	"
" Note:    n1 and nt are acceptable aliases for ns.			"
" 									"
" Example:								"
" suaddhead ns=1024 <bare_traces | sushw key=dt a=4000 >segy_traces	"
" 									"
" This command line adds headers with ns=1024 samples.  The second part	"
" of the pipe sets the trace header field dt to 4 ms.	See also the	"
" selfdocs of related programs  sustrip and supaste.			"
" See:   sudoc supaste							"
" Related Programs:  supaste, sustrip 					"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suaddhead {

/**************** end self doc *******************************************/

/* Credits:
 *	SEP: Einar Kjartansson   c. 1985
 *	CWP: Jack K. Cohen      April 1990
 *      UNOCAL: Zhiming Li	add ascii and binary headers
 */


extern unsigned char su_text_hdr[3200];
extern bhed su_binary_hdr;

segy tr;

void* main_suaddhead( void* args )
{
	int ns;			/* number of samples			*/
	int ftn;		/* ftn=1 for Fortran			*/
	char junk[ISIZE];	/* to discard ftn junk  		*/
	cwp_Bool isreading=cwp_true;    /* true/false flag for while    */

	int ihead=0;		/* counter */
	int iread=0;		/* counter */
	int tsort, ntrpr;	/* Unocal header fields */
#ifdef SU_LINE_HEADER
	cwp_String head;	/* name of file holding headers         */
	FILE *infp=stdin, *outfp=stdout; /* input and output files 	*/
#endif
	FILE *headfp=NULL;	/* . file pointer for pointers		*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suaddhead );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get parameters */
	if (!parObj.getparint("n1", &ns)
	 && !parObj.getparint("nt", &ns)
	 && !parObj.getparint("ns", &ns))  throw cseis_geolib::csException("must specify ns=");
	if (!parObj.getparint("ftn", &ftn))	ftn = 0;
	if (ftn != 0 && ftn != 1)	throw cseis_geolib::csException("ftn=%d must be 0 or 1", ftn);
	if (!parObj.getparint("ntrpr", &ntrpr)) ntrpr = 1;
	if (!parObj.getparint("tsort", &tsort)) tsort = 3;

#ifdef SU_LINE_HEADER

	if (!parObj.getparstring("head"  , &head)) {
		ihead = 0;

	} else {
		ihead = 1;
		if( !(headfp=efopen(head, "r")) ){

                   throw cseis_geolib::csException( "unable to open header file " );
                }
              
	}
        parObj.checkpars();

	/* create id headers */
	if(ihead==0) {
		su_binary_hdr.format = 1;
		su_binary_hdr.ntrpr = ntrpr;
		su_binary_hdr.tsort = tsort;
		su_binary_hdr.fold = ntrpr;
	} else {
		fgethdr(headfp,&su_text_hdr,&su_binary_hdr);
	}

	su_binary_hdr.hns = ns;

#endif
		
	while (isreading==cwp_true) {
		static int tracl = 0;	/* one-based trace number */

		/* If Fortran data, read past the record size bytes */
		if (ftn) efread(junk, ISIZE, 1, stdin);

		/* Do read of data for the segy */
		iread = fread((char *) tr.data, FSIZE, ns, stdin);
		if(iread!=ns) {
			su2cs->setEOF();
			pthread_exit(NULL);
			return retPtr;

		} else {
			if(ihead==0) {
				tr.tracl = ++tracl;
			} else {
				efread(&tr, 1, HDRBYTES, headfp);
			}
			tr.ns = ns;
			tr.trid = TREAL;
			su2cs->putTrace(&tr);
		}

		/* If Fortran data, read past the record size bytes */
		if (ftn) efread(junk, ISIZE, 1, stdin);
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
