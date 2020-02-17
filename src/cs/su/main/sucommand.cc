/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* Copyright (c) Virginia Tech */
/* All rights reserved.		       */

/* SUCOMMAND: $Revision: 1.6 $ ; $Tue Aug 31 11:16:11 EDT 1999 $	*/

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
std::string sdoc_sucommand =
" 									"
" SUCOMMAND - pipe traces having the same key header word to command	"
" 									"
"     sucommand <stdin >stdout [Optional parameters]			"
" 									"
" Required parameters:							"
" 	none								"
" 									"
" Optional parameters: 							"
" 	verbose=0		wordy output				"
" 	key=cdp			header key word to pipe on		"
" 	command=\"suop nop\"    command piped into			"
"	dir=0		0:  change of header key			"
"			-1: break in monotonic decrease of header key	"
"			+1: break in monotonic increase of header key	"
" 									"
" 									"
"Notes:									"
" This program permits limited parallel processing capability by opening"
" pipes for processes, signaled by a change in key header word value.	"
" 									"
" 									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sucommand {


/* Credits:
 *	VT: Matthias Imhof
 *
 * Note:
 *	The "valxxx" subroutines are in su/lib/valpkge.c.  In particular,
 *      "valcmp" shares the annoying attribute of "strcmp" that
 *		if (valcmp(type, val, valnew) {
 *			...
 *		}
 *	will be performed when val and valnew are different.
 *
 */
/**************** end self doc ***********************************/


segy trace;

void* main_sucommand( void* args )
{
	cwp_String key;	/* header key word from segy.h		*/
	cwp_String command;	/* pipe command			*/
	cwp_String type;/* ... its type				*/
	int indx;	/* ... its index			*/
	int nsegy;	/* number of bytes in the segy		*/
	Value val;	/* value of key in current gather	*/
	Value valnew;	/* value of key in trace being treated	*/
	int fold;	/* number of traces per gather		*/
	int verbose;	/* verbose flag				*/
	int dir;	/* direction flag			*/
	FILE *pipe;

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sucommand );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Set parameters */
	if (!parObj.getparint   ("verbose", &verbose))	 verbose = 0;
	if (!parObj.getparstring("command", &command))	 command = "suop nop";
	if (!parObj.getparstring("key", &key))		 key = "cdp";
	if (!parObj.getparint   ("dir", &dir))	 	 dir = 0;
        parObj.checkpars();
	dir=-dir;

	type = hdtype(key);
	indx = getindex(key);

	/* Set up for first trace (must compare new key field each time) */
	pipe = (FILE *) popen(command, "w");
	fold = 1;
	
	nsegy = cs2su->getTrace(&trace);
	gethval(&trace, indx, &val);
	fputtr(pipe, &trace);

	/* Loop over traces */
	for(;;)
	{	int	cmp;
		nsegy = cs2su->getTrace(&trace);
		if (! nsegy) break;
		
		gethval(&trace, indx, &valnew);
		cmp = valcmp(type, val, valnew);
		val = valnew;
		
		if (((dir==0) && cmp) || ((dir != 0) && (dir!=cmp))) {	
			/* val and valnew differ, indicating a  */
			/* new gather ? 				*/
			if (verbose) {
				fprintf(stderr, "val=");
				fprintfval(stderr, type, val);
				fprintf(stderr, "\tgather size=%d\n", fold);
			}

			pclose(pipe);

			/* Set up for next gather */
			
			pipe = (FILE *) popen(command, "w");
			fold = 1;
			fputtr(pipe, &trace);
		}
		
		else {	/* still in same gather */
			fputtr(pipe, &trace);
			++fold;
		}
	}
	
	if (verbose) {  fprintf(stderr, "val=");
			fprintfval(stderr, type, val);
			fprintf(stderr, "\tgather size=%d\n", fold);
	}
	pclose(pipe);
	
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
