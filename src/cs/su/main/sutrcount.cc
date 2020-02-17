/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUTRCOUNT: $Revision: 1.4 $ ; $Date: 2012/01/11 00:44:40 $	*/


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
extern "C" {
  #include <stdio.h>

}
extern "C" {
  #include <stdlib.h>

}
extern "C" {
  #include <math.h>

}
#include "su.h"
#include "segy.h"

/*********************** self documentation *****************************/
std::string sdoc_sutrcount =
" SUTRCOUNT - SU program to count the TRaces in infile		"
"       							"
"   sutrcount < infile					     	"
" Required parameters:						"
"       none							"
" Optional parameter:						"
"    outpar=stdout						"
" Notes:       							"
" Once you have the value of ntr, you may set the ntr header field"
" via:      							"
"       sushw key=ntr a=NTR < datain.su  > dataout.su 		"
" Where NTR is the value of the count obtained with sutrcount 	"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sutrcount {


/*
 * Credits:  B.Nemeth, Potash Corporation, Saskatchewan 
 * 		given to CWP in 2008 with permission of Potash Corporation
 */

/**************** end self doc ********************************/
   
/* Segy data constants */
segy tr;				/* SEGY trace */

void* main_sutrcount( void* args )
{
	/* Segy data constans */
	int ntr=0;		/* number of traces			*/
	char *outpar=NULL;	/* name of file holding output		*/
	FILE *outparfp=stdout;	/* ... its file pointer			*/

	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sutrcount );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	
	/* Get information from the first header */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	if (!parObj.getparstring("outpar", &outpar))	outpar = "/dev/stdout" ;
	
	outparfp = efopen(outpar, "w");

        parObj.checkpars();
	/* Loop over traces getting a count */
	do {
		++ntr;
	} while(cs2su->getTrace(&tr));

	fprintf(outparfp, "%d", ntr);

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
