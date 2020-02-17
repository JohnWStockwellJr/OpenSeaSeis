/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSPIKE: $Revision: 1.15 $ ; $Date: 2011/11/12 00:40:42 $	*/

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
std::string sdoc_suspike =
"								"
" SUSPIKE - make a small spike data set 			"
"								"
" suspike [optional parameters] > out_data_file  		"
"								"
" Creates a common offset su data file with up to four spikes	"
" for impulse response studies					"
"								"
" Optional parameters:						"
"	nt=64 		number of time samples			"
"	ntr=32		number of traces			"
" 	dt=0.004 	time sample rate in seconds		"
" 	offset=400 	offset					"
"	nspk=4		number of spikes			"
"	ix1= ntr/4	trace number (from left) for spike #1	"
"	it1= nt/4 	time sample to spike #1			"
"	ix2 = ntr/4	trace for spike #2			"
"	it2 = 3*nt/4 	time for spike #2			"
"	ix3 = 3*ntr/4;	trace for spike #3			"
"	it3 = nt/4;	time for spike #3			"
"	ix4 = 3*ntr/4;	trace for spike #4			"
"	it4 = 3*nt/4;	time for spike #4			"
"								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suspike {


/* Credits:
 *	CWP: Shuki Ronen, Chris Liner
 *
 * Trace header fields set: ns, dt, offset
 */
/**************** end self doc ***********************************/


segy tr;

void* main_suspike( void* args )
{
	int nt;				/* number of time samples	*/
	int ntr;			/* number of traces		*/
	int itr;			/* trace counter		*/
	int nspk;			/* number of spikes		*/
	int it1;			/* time of 1st spike		*/
	int ix1;			/* position of 1st spike	*/
	int it2;			/* time of 2nd spike		*/
	int ix2;			/* position of 2nd spike	*/
	int ix3;			/* position of 3rd spike	*/
	int it3;			/* time of 3rd spike		*/
	int ix4;			/* position of 4th spike	*/
	int it4;			/* time of 4th spike		*/
	float dt;			/* time sampling interval	*/
	float offset;			/* offset			*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suspike );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	nt = 64;	parObj.getparint("nt", &nt);
	CHECK_NT("nt",nt);				tr.ns = nt;
	ntr = 32;	parObj.getparint("ntr", &ntr);
	dt = 0.004;	parObj.getparfloat("dt", &dt);		tr.dt = dt*1000000;
	offset = 400;	parObj.getparfloat("offset", &offset);	tr.offset = offset;
	nspk = 4;	parObj.getparint("nspk", &nspk);
	ix1 = ntr/4;	parObj.getparint("ix1", &ix1); 
	it1 = nt/4;	parObj.getparint("it1", &it1);
	ix2 = ntr/4;	parObj.getparint("ix2", &ix2);
	it2 = 3*nt/4;	parObj.getparint("it2", &it2);
	ix3 = 3*ntr/4;	parObj.getparint("ix3", &ix3);
	it3 = nt/4;	parObj.getparint("it3", &it3);
	ix4 = 3*ntr/4;	parObj.getparint("ix4", &ix4);
	it4 = 3*nt/4;	parObj.getparint("it4", &it4);

	for (itr = 0; itr < ntr; itr++) {
		memset( (void *) tr.data, 0, nt * FSIZE);
		if (itr == ix1-1) tr.data[it1-1] = 1.0;  
		if (nspk > 1 && itr == ix2-1) tr.data[it2-1] = 1.0;
		if (nspk > 2 && itr == ix3-1) tr.data[it3-1] = 1.0;
		if (nspk > 3 && itr == ix4-1) tr.data[it4-1] = 1.0;
		tr.tracl = itr + 1;
		su2cs->putTrace(&tr);
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
