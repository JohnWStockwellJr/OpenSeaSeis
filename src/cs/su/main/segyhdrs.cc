/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SEGYHDRS: $Revision: 1.31 $ ; $Date: 2011/11/16 17:43:20 $	*/

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
#include "su_xdr.h"
#include "bheader.h"
#include "tapesegy.h"

/*********************** self documentation **********************/
std::string sdoc_segyhdrs =
" 									"
" SEGYHDRS - make SEG-Y ascii and binary headers for segywrite		"
" 									"
" segyhdrs [ < sudata ] [optional parameters] [ > copy of sudata ]      "
" 									"
" Required parameters:							"
"	ns=  if no input trace header					"
"	dt=  if no input trace header					"
" Optional parameters:							"
" 	ns=tr.ns from header    number of samples on input traces	"
" 	dt=tr.dt from header	sample rate (microseconds) from traces	"
" 	bfile=binary		name of file containing binary block	"
" 	hfile=header		name of file containing ascii block	"
"   Some binary header fields are set:					"
" 	jobid=1			job id field				"
" 	lino=1			line number (only one line per reel)	"
" 	reno=1			reel number				"
" 	format=1		data format				"
" 									"
" All other fields are set to 0, by default.				"
" To set any binary header field, use sukeyword to find out		"
" the appropriate keyword, then use the getpar form:			"
" 	keyword=value	to set keyword to value				"
" 									"
" The header file is created as ascii and is translated to ebcdic	"
" by segywrite before being written to tape.  Its contents are		"
" formal but can be edited after creation as long as the forty		"
" line format is maintained.						"
" 									"
" Caveat: This program has not been tested under XDR for machines       "
"	 not having a 2 byte unsigned short integral data type.	"
" 									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace segyhdrs {


/* Credits:
 *
 *	CWP: Jack K. Cohen,  John Stockwell 
 *      MOBIL: Stew Levin
 */
/**************** end self doc ***********************************/

static segy tr;
static bhed bh;

void* main_segyhdrs( void* args )
{
#ifdef SUXDR
	XDR  bhed_xdr;		/* xdr structure to associate w/binaryfp */
#endif
	char *hfile;		/* name of ascii header file		*/
	char buf[81];		/* buffer for ascii file lines + null	*/
	char *bfile;		/* name of binary header file		*/
	FILE *headerfp;		/* file pointer for hfile		*/
	FILE *binaryfp;		/* file pointer for bfile		*/
	int i;			/* counter				*/
	cwp_Bool redin=cwp_false;
	cwp_Bool redout=cwp_false; /* flags for redirected input/output    */

	/* bhed fields the user might wish to set, taken right out of segy.h */
	int jobid;	/* job identification number */
	int lino;	/* line number (only one line per reel) */
	int reno;	/* reel number */
	short ntrpr;	/* number of data traces per record */
	short nart;	/* number of auxiliary traces per record */
	unsigned short hdt; /* sample interval in microsecs for this reel */
	unsigned short dto; /* same for original field recording */

	/* hns replaces ns in original version of segyhdrs */
	unsigned short hns; /* number of samples per trace for this reel */

	unsigned short nso; /* same for original field recording */
	short format;	/* data sample format code:
				1 = floating point (4 bytes)
				2 = fixed point (4 bytes)
				3 = fixed point (2 bytes)
				4 = fixed point w/gain code (4 bytes) */
	short fold;	/* CDP fold expected per CDP ensemble */
	short tsort;	/* trace sorting code: 
				1 = as recorded (no sorting)
				2 = CDP ensemble
				3 = single fold continuous profile
				4 = horizontally stacked */
	short vscode;	/* vertical sum code:
				1 = no sum
				2 = two sum ...
				N = N sum (N = 32,767) */
	short hsfs;	/* sweep frequency at start */
	short hsfe;	/* sweep frequency at end */
	short hslen;	/* sweep length (ms) */
	short hstyp;	/* sweep type code:
				1 = linear
				2 = parabolic
				3 = exponential
				4 = other */
	short schn;	/* trace number of sweep channel */
	short hstas;	/* sweep trace taper length at start if
			   tapered (the taper starts at zero time
			   and is effective for this length) */
	short hstae;	/* sweep trace taper length at end (the ending
			   taper starts at sweep length minus the taper
			   length at end) */
	short htatyp;	/* sweep trace taper type code:
				1 = linear
				2 = cos-squared
				3 = other */
	short hcorr;	/* correlated data traces code:
				1 = no
				2 = yes */
	short bgrcv;	/* binary gain recovered code:
				1 = yes
				2 = no */
	short rcvm;	/* amplitude recovery method code:
				1 = none
				2 = spherical divergence
				3 = AGC
				4 = other */
	short mfeet;	/* measurement system code:
				1 = meters
				2 = feet */
	short polyt;	/* impulse signal polarity code:
				1 = increase in pressure or upward
				    geophone case movement gives
				    negative number on tape
				2 = increase in pressure or upward
				    geophone case movement gives
				    positive number on tape */
	short vpol;	/* vibratory polarity code:
				code	seismic signal lags pilot by
				1	337.5 to  22.5 degrees
				2	 22.5 to  67.5 degrees
				3	 67.5 to 112.5 degrees
				4	112.5 to 157.5 degrees
				5	157.5 to 202.5 degrees
				6	202.5 to 247.5 degrees
				7	247.5 to 292.5 degrees
				8	293.5 to 337.5 degrees */


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_segyhdrs );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	redin = cwp_true;
	switch(filestat(STDIN)) {
	case TTY:
	case DIRECTORY:
	case BADFILETYPE:
		redin = cwp_false;
	break;
	default: /* Others OK */
	break;
	}
	redout = cwp_true;
	switch(filestat(STDOUT)) {
	case TTY:
	case DIRECTORY:
	case BADFILETYPE:
		redout = cwp_false;
	break;
	default: /* Others OK */
	break;
	}

	if(redin==cwp_false)

	
	/* get number of samples and sample rate from first trace or parObj.getpar */
	if (redin==cwp_true) {
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	if (!parObj.getparushort("ns", &hns))		hns = tr.ns;
	if (!parObj.getparushort("dt", &hdt))		hdt = tr.dt;
	if (!hns) warn("hns not set in binary header, consider using ns=");
	if (!hdt) warn("hdt not set in binary header, consider using dt=");
	} else {
	    if(!parObj.getparushort("ns",&hns)) throw cseis_geolib::csException("need ns=");
	    if(!parObj.getparushort("dt",&hdt)) throw cseis_geolib::csException("need dt=");
	}

	/* Get parameters */
	if (!parObj.getparstring("hfile", &hfile))	hfile = "header";
	if (!parObj.getparstring("bfile", &bfile))	bfile = "binary";

	/* Optional settings; user has access to all named fields of bhed */
	if (!parObj.getparint("jobid", &jobid))	jobid = 1;
	if (!parObj.getparint("lino", &lino))		lino = 1;
	if (!parObj.getparint("reno", &reno))		reno = 1;
	if (!parObj.getparshort("ntrpr", &ntrpr))	ntrpr = 0;
	if (!parObj.getparshort("nart", &nart))	nart = 0;
	if (!parObj.getparushort("dto", &dto))		dto = 0;
	if (!parObj.getparushort("nso", &nso))		nso = 0;
	if (!parObj.getparshort("format", &format))	format = 1;
	if (!parObj.getparshort("fold", &fold))	fold = 0;
	if (!parObj.getparshort("tsort", &tsort))	tsort = 0;
	if (!parObj.getparshort("vscode", &vscode))	vscode = 0;
	if (!parObj.getparshort("hsfs", &hsfs))	hsfs = 0;
	if (!parObj.getparshort("hsfe", &hsfe))	hsfe = 0;
	if (!parObj.getparshort("hslen", &hslen))	hslen = 0;
	if (!parObj.getparshort("hstyp", &hstyp))	hstyp = 0;
	if (!parObj.getparshort("schn", &schn))	schn = 0;
	if (!parObj.getparshort("hstas", &hstas))	hstas = 0;
	if (!parObj.getparshort("hstae", &hstae))	hstae = 0;
	if (!parObj.getparshort("htatyp", &htatyp))	htatyp = 0;
	if (!parObj.getparshort("hcorr", &hcorr))	hcorr = 0;
	if (!parObj.getparshort("bgrcv", &bgrcv))	bgrcv = 0;
	if (!parObj.getparshort("rcvm", &rcvm))	rcvm = 0;
	if (!parObj.getparshort("mfeet", &mfeet))	mfeet = 0;
	if (!parObj.getparshort("polyt", &polyt))	polyt = 0;
	if (!parObj.getparshort("vpol", &vpol))	vpol = 0;
        parObj.checkpars();


	/* Open files for writing */
	headerfp = efopen(hfile, "w");
	binaryfp = efopen(bfile, "w");
#ifdef SUXDR
	xdrstdio_create(&bhed_xdr, binaryfp, XDR_ENCODE);
#endif

	/* Create ascii header */
	sprintf(buf, "%-79.79s\n", "C      This tape was made at the");
	efwrite((char *) buf, 1, 80, headerfp);

	sprintf(buf, "%-79.79s\n", "C");
	efwrite((char *) buf, 1, 80, headerfp);

	sprintf(buf, "%-79.79s\n", "C      Center for Wave Phenomena");
	efwrite((char *) buf, 1, 80, headerfp);

	sprintf(buf, "%-79.79s\n", "C      Colorado School of Mines");
	efwrite((char *) buf, 1, 80, headerfp);

	sprintf(buf, "%-79.79s\n", "C      Golden, CO, 80401");
	efwrite((char *) buf, 1, 80, headerfp);

	for (i = 0;  i < 35; i++) {
		sprintf(buf, "%-79.79s\n", "C");
		efwrite((char *) buf, 1, 80, headerfp);
	}


	/* Create binary header; set all named fields */
	memset((void *) &bh, 0, BNYBYTES);
	bh.jobid = (int) jobid;
	bh.lino = (int) lino;
	bh.reno = (int) reno;
	bh.ntrpr = (short) ntrpr;
	bh.nart = (short) nart;
	bh.hdt = (short) hdt;
	bh.dto = (short) dto;
	bh.hns = (short) hns;
	bh.nso = (short) nso;
	bh.format = (short) format;
	bh.fold = (short) fold;
	bh.tsort = (short) tsort;
	bh.vscode = (short) vscode;
	bh.hsfs = (short) hsfs;
	bh.hsfe = (short) hsfe;
	bh.hslen = (short) hslen;
	bh.hstyp = (short) hstyp;
	bh.schn = (short) schn;
	bh.hstas = (short) hstas;
	bh.hstae = (short) hstae;
	bh.htatyp = (short) htatyp;
	bh.hcorr = (short) hcorr;
	bh.bgrcv = (short) bgrcv;
	bh.rcvm = (short) rcvm;
	bh.mfeet = (short) mfeet;
	bh.polyt = (short) polyt;
	bh.vpol = (short) vpol;

	/* Write binary header from bh structure to designated file */
#ifdef SUXDR
	if(FALSE == xdrbhdrsub(&bhed_xdr,&bh))
		throw cseis_geolib::csException("%s: problem writing binary tape header",__FILE__);
#else
	efwrite( (char *) &bh, 1, BNYBYTES, binaryfp);
#endif

	/* Clean up */
	efclose(headerfp);
#ifdef SUXDR
	xdr_destroy(&bhed_xdr);
#endif
	efclose(binaryfp);

	if(redin==cwp_true) {
		/* suck in all input traces to avoid broken pipe */
	  do {
	      if(redout==cwp_true) su2cs->putTrace(&tr);
	  } while(cs2su->getTrace(&tr));
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
