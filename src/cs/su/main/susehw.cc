/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSETEHW: $Revision: 1.3 $ ; $Date: 2012/03/28 16:44:34 $	*/

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
std::string sdoc_susehw =
"									"
" SUSEHW - Set the value the Header Word denoting trace number within	"
"	     an Ensemble defined by the value of another header word	"
"									"
"     susehw <stdin >stdout [options]					"
"									"
" Required Parameters:							"
"	none								"
"									"
" Optional Parameters:							"
" key1=cdp	Key header word defining the ensemble			"
" key2=cdpt	Key header word defining the count within the ensemble	"
" a=1		starting value of the count in the ensemble		"
" b=1		increment or decrement within the ensemble		"
"									"
" Notes:								"
" This code was written because suresstat requires cdpt to be set.	"
" The computation is 							"
" 	val(key2) = a + b*i						"
"									"
" The input data must first be sorted into constant key1 gathers.	"
" Example: setting the cdpt field					" 
"        susetehw < cdpgathers.su a=1 b=1 key1=cdp key2=cdpt > new.su	"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace susehw {


/* Credits:
 *  CWP: John Stockwell (Feb 2008) in answer to a question by Warren Franz
 *        based on various codes, including susplit, susshw, suchw
 */
/**************** end self doc *******************************************/

/* Prototypes of functions used internally */
void setval(cwp_String type, Value *valp, double a, double b, double i);

segy tr;

void* main_susehw( void* args )
{
	cwp_String key1;	/* key word labeling the ensemble	*/
	Value val1;		/* value of key1			*/
	int ival1;		/* key1 value as integer 		*/
	int ival1old;		/* last value of ival			*/ 
	cwp_String type1;	/* key1's type				*/
	int index1;		/* index of key1			*/

	int i=0;		/* counter 				*/

	cwp_String key2;	/* key word labelin ensemble members	*/
	Value val2;		/* value of key2			*/
	cwp_String type2;	/* key2's type				*/
	int index2;		/* index of key2			*/

	double a;		/* initial header value within ensemble */
	double b;		/* increment (decrement) within ensemble*/

	int verbose;		/* =1 to echo filenames, etc.		*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_susehw );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Default parameters;  User-defined overrides */
	if (!parObj.getparstring("key1", &key1))	key1="cdp";
	if (!parObj.getparstring("key2", &key2))	key2="cdpt";
	if (!parObj.getpardouble("a", &a))	 	a=1;
	if (!parObj.getpardouble("b", &b))	 	b=1;
	if (!parObj.getparint("verbose",&verbose)) 	verbose=0;

        parObj.checkpars();

	/* Evaluate time bounds from parObj.getpars and first header */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");

	/* Set up key values */
	type1 = hdtype(key1);
	index1 = getindex(key1);
	type2 = hdtype(key2);
	index2 = getindex(key2);

	/* Output debugging information */
	if (verbose) warn("key1 = %s",key1);
	if (verbose) warn("key2 = %s",key2);

	/* Get value of key1 from the first trace */
	gethval(&tr, index1, &val1);
	ival1 = vtoi(type1, val1);
	ival1old = ival1;

	/* Main loop over traces */
	i = 0;
	do {
		/* get value of key1 header field */
		gethval(&tr, index1, &val1);
		ival1 = vtoi(type1, val1);

		/* test to see if we are in a new ensemble */
		if (ival1!=ival1old) { 
			ival1old = ival1;
			i = 0; /* beginning of a new ensemble */
		}

		/* set header value in trace */
		setval(type2, &val2, a,  b,  i);
		puthval(&tr,index2,&val2);

		/* increment trace counter */
		++i;

		/* output trace with modified header */
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

void setval( cwp_String type, Value *valp, double a, double b, double i)
/**********************************************************************
sethval - set header word value on a trace
***********************************************************************
cwp_String type     type (int, short, long, double, float)
Value *valp         pointer to value
double a            initial header value
double b            increment (decrement) of header values
double i	    trace counter
***********************************************************************
Notes:
The SU header fields have a variety of types. Getting this right requires
a bit of doing.  
***********************************************************************
Author:  CWP: John Stockwell,  Feb 2008, hacked from the one in sushw.
**********************************************************************/
{
	switch (*type) {
	case 's':
		throw cseis_geolib::csException("can't set char header word");
	break;
	case 'h':
		valp->h = (a + b * i);
	break;
	case 'u':
		valp->u = (a + b * i);
	break;
	case 'l':
		valp->l = (long) (a + b * i );
	break;
	case 'v':
		valp->v = (unsigned long) (a + b * i);
	break;
	case 'i':
		valp->i = (a + b * i );
	break;
	case 'p':
		valp->p = (a + b * i );
	break;
	case 'f':
		valp->f = (a + b * i );
	break;
	case 'd' :
		valp->d = (a + b * i );
	default:
		throw cseis_geolib::csException("unknown type %s", type);
	break;
	}
	return;
}


} // END namespace
