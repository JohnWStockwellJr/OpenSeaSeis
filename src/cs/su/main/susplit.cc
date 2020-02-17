/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSPLIT */

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
#include "bheader.h"
#include "tapebhdr.h"


/*********************** self documentation ******************************/
std::string sdoc_susplit =
"									"
" SUSPLIT - Split traces into different output files by keyword value	"
"									"
"     susplit <stdin >stdout [options]					"
"									"
" Required Parameters:							"
"	none								"
"									"
" Optional Parameters:							"
"	key=cdp		Key header word to split on (see segy.h)	"
"	stem=split_	Stem name for output files			"
"	middle=key	middle of name of output files			"
"	suffix=.su	Suffix for output files				"
"	numlength=7	Length of numeric part of filename		"
"	verbose=0	=1 to echo filenames, etc.			"
"	close=1		=1 to close files before opening new ones	"
"									"
" Notes:								"
" The most efficient way to use this program is to presort the input data"
" into common keyword gathers, prior to using susplit.			"
"									"
" Use \"suputgthr\" to put SU data into SU data directory format.	"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace susplit {


/* Credits:
 *	Geocon: Garry Perratt hacked together from various other codes
 * 
 */
/**************** end self doc *******************************************/


segy tr;

void* main_susplit( void* args )
{
	cwp_String key;		/* header key word from segy.h		*/
	int ns;			/* ns as an  int			*/
	int numlength;		/* length of split key number format	*/
	char format[BUFSIZ];	/* output filename format		*/

	int index;		/* index of key				*/
	Value val;		/* value of key				*/
	int val_i=0;		/* key value as integer 		*/
	int lastval_i=0;	/* last key value as integer	 	*/

	cwp_String type;	/* key's type				*/
	char filename[BUFSIZ];	/* output file name			*/

	cwp_String stem;	/* output file stem			*/
	cwp_String middle;	/* output file middle			*/
	cwp_String suffix;	/* output file suffix			*/
	FILE *outfp=NULL;	/* pointer to output file		*/

	int init;		/* initialisation flag for first efopen	*/
	int verbose;		/* =1 to echo filenames, etc.		*/
	int close;	/* =1 to close files before opening a new one	*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_susplit );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	init=0;

	/* Default parameters;  User-defined overrides */
	if (!parObj.getparstring("key", &key))	 key="cdp";
	if (!parObj.getparstring("stem", &stem)) stem="split_";
	if (!parObj.getparstring("middle", &middle)) middle=key;
	if (!parObj.getparstring("suffix", &suffix)) suffix=".su";
	if (!parObj.getparint("numlength", &numlength)) numlength=7;
	if (!parObj.getparint("verbose",&verbose)) verbose=0;
	if (!parObj.getparint("close",&close)) close=1;
        parObj.checkpars();

	/* Evaluate time bounds from parObj.getpars and first header */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	ns = tr.ns;

	type = hdtype(key);
	index = getindex(key);

	if (verbose) warn("key is %s",key);

	/* Main loop over traces */
	do {
		gethval(&tr, index, &val);
		val_i = vtoi(type, val);
		
		if (val_i!=lastval_i || init==0) {
			if (init==0) init=1;
			else if (close) efclose(outfp);
			(void)sprintf(format, "%%s%%s%%0%dd%%s",numlength);
			(void)sprintf(filename, format, stem, middle, val_i, suffix);
			if (verbose) warn("output file is %s",filename);
			outfp = efopen(filename, "ab");
		}
		fputtr(outfp,&tr);
		lastval_i = val_i;

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
