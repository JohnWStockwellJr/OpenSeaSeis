/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUTAPER: $Revision: 1.5 $ ; $Date: 2011/11/16 23:33:10 $	*/

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
std::string sdoc_sugausstaper =
" 								"
" SUGAUSSTAPER - Multiply traces with gaussian taper		"
" 								"
" sugausstaper < stdin > stdout [optional parameters]		"
" 								"
" Required Parameters:					   	"
"   <none>							"
"								"
" Optional parameters:					   	"
" key=offset    keyword of header field to weight traces by 	"
" x0=300        key value defining the center of gaussian window" 
" xw=50         width of gaussian window in units of key value 	"
"								"
" Notes:							"
" Traces are multiplied with a symmetrical gaussian taper 	"
"  	w(t)=exp(-((key-x0)/xw)**2)				"
" unlike \"sutaper\" the value of x0 defines center of the taper"
" rather than the edges of the data.				"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sugausstaper {


/* Credits:
 *
 *	Thomas Bohlen, formerly of TU Bergakademie, Freiberg GDR
 *      most recently of U Karlsruhe
 *          04.01.2002
 *
 * Trace header fields accessed: ns
 */
/**************** end self doc ***********************************/


segy tr;

void* main_sugausstaper( void* args )
{
	char *key;	/* header key word from segy.h		*/
	char *type;     /* ... its type				*/	
	int index;	/* ... its index			*/
	Value val;	/* ... its value			*/
	float fval;     /* ... its value cast to float		*/

	int ns;		/* number of sample points on traces	*/
	float x0, xw;	/* centre and width of gauss taper	*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sugausstaper );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get optional parameters */
	if (!parObj.getparstring("key", &key))	key = "offset";
	if (!parObj.getparfloat("x0",&x0))	x0 = 300.0;
	if (!parObj.getparfloat("xw",&xw))	xw = 50.0;

        parObj.checkpars();
	/* Get key type and index */
	type = hdtype(key);
	index = getindex(key);

	/* Get info from first trace */
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	ns = tr.ns;
 
	/* Loop through traces */
	do {
		register int i = 0;	     /* counter */

		/* Get value of key and convert to float */
		gethval(&tr, index, &val);
		fval = vtof(type,val);

		/* Loop over samples in trace and apply weighting */
		for (i=0; i < ns; ++i)
			tr.data[i] *= exp(-pow((fval-x0)/xw,2.0));

		/* Put out weighted traces */
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
