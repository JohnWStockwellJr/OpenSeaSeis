/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* suswapbytes: $Revision: 1.25 $ ; $Date: 2011/11/16 17:43:20 $	*/

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
std::string sdoc_suswapbytes =
"									"
" SUSWAPBYTES - SWAP the BYTES in SU data to convert data from big endian"
"               to little endian byte order, and vice versa		"
"									"
" suswapbytes < stdin [optional parameter] > sdtout			"
"									"
" 	format=0		foreign to native			"
" 				=1 native to foreign			"
" 	ns=from header		if ns not set in header, must be set by hand"
" Notes:								"
"  The 'native'	endian is the endian (byte order) of the machine you are"
"  running this program on. The 'foreign' endian is the opposite byte order."
"   									"
" Examples of big endian machines are: IBM RS6000, SUN, NeXT		"
" Examples of little endian machines are: PCs, DEC			"
"									"
" Caveat: this code has not been tested on DEC				"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suswapbytes {


/* Credits: 
 *	CWP: adapted for SU by John Stockwell 
 *		based on a code supplied by:
 *	Institute fur Geophysik, Hamburg: Jens Hartmann (June 1993)
 *
 * Trace header fields accessed: ns
 */
/**************** end self doc ********************************/

/* function prototypes for subroutines used internally */

static segy tr;

void* main_suswapbytes( void* args )
{
	size_t trbytes;
	size_t databytes;
	int format,ns,i;
	float *data;

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suswapbytes );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Make sure stdout is a file or pipe */
	switch(filestat(STDOUT)) {
	case TTY:
		throw cseis_geolib::csException("stdout can't be tty");
	break;
	case DIRECTORY:
		throw cseis_geolib::csException("stdout must be a file, not a directory");
	break;
	case BADFILETYPE:
		throw cseis_geolib::csException("stdout is illegal filetype");
	break;
	default:
			/* all others ok */
	break;
	}


	if (!parObj.getparint("format", &format)) format = 0;

	if (format==0) { /* convert foreign to native byte order */

		/* read the header off the first trace */
		efread(&tr,1,HDRBYTES,stdin);

		/* swap the header values on the first trace */
		for (i = 0; i < SU_NKEYS; ++i) swaphval(&tr,i);

		if (!parObj.getparint("ns", &ns))	ns = tr.ns ;
                parObj.checkpars();
 
		if (!ns) throw cseis_geolib::csException("ns not set in header, set ns by hand!");

		/* size of data portion */
		databytes = ns*FSIZE;
		
		/* allocate space for data */
		data = alloc1float(databytes);

		/* read the data values on the first trace */
		efread(data,1,databytes,stdin);

		/* now swap the data values */
		for (i = 0; i < ns ; ++i) swap_float_4(&data[i]);

		/* write out the first trace */
		efwrite(&tr, 1, HDRBYTES, stdout);
		efwrite(data, 1, databytes, stdout);

		/* total size of a trace */
		trbytes= HDRBYTES + databytes;

		/* read, swap, and output the rest of the traces */
 		while (efread(&tr,1,trbytes,stdin)) {

			/* swap header values */
			for (i = 0; i < SU_NKEYS; ++i) swaphval(&tr,i);

			/* now swap the data values */
			for (i = 0; i < ns ; ++i)
					swap_float_4(&tr.data[i]);

			efwrite(&tr, 1, trbytes, stdout);
		} /* finished converting foreign to native */

	} else { /* convert "native" to "foreign" byte order */
                parObj.checkpars();
		/* read in traces, swap, and write to stdout */
		while (cs2su->getTrace(&tr)) {

			trbytes=HDRBYTES+FSIZE*tr.ns;

			/* now swap the data values */
			for (i = 0; i < ((int) tr.ns) ; ++i)
					swap_float_4(&tr.data[i]);

			/* swap header values */
			for (i = 0; i < SU_NKEYS; ++i) swaphval(&tr,i);


			efwrite(&tr,1,trbytes,stdout);
		} /* finished converting native to foreign */
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
