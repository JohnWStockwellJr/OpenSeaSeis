/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUFNZERO: $Revision: 1.5 $ ; $Date: 2011/12/23 19:36:06 $	*/

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

/*************************** self documentation **************************/
std::string sdoc_sufnzero =
"									"
" SUFNZERO - get Time of First Non-ZERO sample by trace              "
"									"
"  sufnzero <stdin >stdout [optional parameters] 			"
"									"
" Required parameters:							"
"	none								"
"									"
" Optional parameters: 							"
"	mode=first   	Output time of first non-zero sample		"
"	             	=last for last non-zero sample			"
"	             	=both for both first & last non-zero samples    "
"									"
"	min=0   	Threshold value for considering as zero         "
"			Any abs(sample)<min is considered to be zero	"
"									"
"	key=key1,...	Keyword(s) to print				"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sufnzero {


/* Credits:
 *      Geocon : Garry Perratt
 *	based on surms by the same, itself based on sugain & sumax by:
 *	CWP : John Stockwell
 */
/**************** end self doc ***********************************/

segy tr;

void* main_sufnzero( void* args )
{
	cwp_String key[SU_NKEYS];	/* array of keywords		*/
	int ikey;			/* input key counter		*/
	int nkeys;			/* number of keywords to be gotten*/
	int index;			/* index of header keys		*/
	Value keyval;			/* header key value		*/
	cwp_String keytype;		/* header key type		*/
	float min=0.0;			/* Threshold value considered zero*/
	char *mode;			/* desired output		*/
	int nt;				/* number of time points on trace*/
	int itr=0;			/* trace number	-- bumped at loop top*/
	int i=0;			/* time point number		*/
	float ftime=0.0;		/* time of first non-zero sample*/
	float ltime=0.0;		/* time of last non-zero sample	*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sufnzero );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get optional parameters */
	if (!parObj.getparstring("mode", &mode)) mode="first";
	if (!(	STREQ(mode, "first") || 
		STREQ(mode, "last") ||
		STREQ(mode, "both")
	)) throw cseis_geolib::csException("%s unknown mode", mode); 

	if (!parObj.getparfloat("min", &min)) min=0.0;

	if ((nkeys=parObj.countparval("key"))!=0) parObj.getparstringarray("key",key);

        parObj.checkpars();
	/* Get info from first trace */
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;

	/* Loop through data */
	do {
		itr++ ;

		ftime=-1.0;
		ltime=-1.0;
		for (i = 0; i < nt; ++i) {
			/* fabs for absolute value of float */
			if (fabs(tr.data[i])>min && ftime==-1) 
					ftime=i*tr.dt/1000;

			if (fabs(tr.data[i])>min)
					ltime=i*tr.dt/1000;
		}

		/* Output ASCII */
		for (ikey = 0; ikey < nkeys; ++ikey) {
			index=getindex(key[ikey]);
			gethval(&tr,index,&keyval);
			keytype=hdtype(ikey[key]);
			printfval(keytype, keyval); 
			printf(" ");
		}
		if (STREQ(mode, "first")) printf("%g\n", ftime); 
		else if (STREQ(mode, "last")) printf("%g\n", ltime);
		else if (STREQ(mode, "both")) printf("%g %g\n", ftime, ltime); 

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
