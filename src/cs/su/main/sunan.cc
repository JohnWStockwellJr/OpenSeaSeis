/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUNAN: $Revision: 1.8 $ ; $Date: 2011/12/21 23:19:56 $        */

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
  #include <math.h>

}
#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
std::string sdoc_sunan =
"								"
" SUNAN - remove NaNs & Infs from the input stream		"
"								"
"    sunan < in.su >out.su					"
"								"
" Optional parameters:						"
" verbose=1	echo locations of NaNs or Infs to stderr	"
"	        =0 silent					"
" ...user defined ... 						"
"								"
" value=0.0	NaNs and Inf replacement value			"
" ... and/or....						"
" interp=0	=1 replace NaNs and Infs by interpolating	"
"                   neighboring finite values			"
"								"
" Notes:							"
" A simple program to remove NaNs and Infs from an input stream."
" The program sets NaNs and Infs to \"value\" if interp=0. When	"
" interp=1 NaNs are replaced with the average of neighboring values"
" provided that the neighboring values are finite, otherwise	"
" NaNs and Infs are replaced by \"value\".			"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sunan {


/*
 * Author: Reginald H. Beardsley  2003   rhb@acm.org
 *
 *  A simple program to remove NaNs & Infs from an input stream. They
 *  shouldn't be there, but it can be hard to find the cause and fix
 *  the problem if you can't look at the data.
 *
 *  Interpolation idea comes from a version of sunan modified by
 *  Balasz Nemeth while at Potash Corporation in Saskatchewan.
 *
 */

/**************** end self doc ********************************/

segy tr;

void* main_sunan( void* args )
{

	int i;			/* counter			*/
	int itr=0;		/* trace counter		*/
	int verbose;		/* =0 silent,  =1 chatty	*/
	int interp;		/* =1 interpolate to get NaN	*/
				/* and Inf replacement values	*/
			
	float value;		/* value to set NaN and Infs to */

	/* Initialize */
   	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
   	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
   	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
   	int argc = suArgs->argc;
   	char **argv = suArgs->argv;
   	cseis_su::csSUGetPars parObj;

   	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
   	su2cs->setSUDoc( sdoc_sunan );
   	if( su2cs->isDocRequestOnly() ) return retPtr;
   	parObj.initargs(argc, argv);

   	try {  /* Try-catch block encompassing the main function body */


	/* Get info from first trace */
	if(!cs2su->getTrace(&tr) ) throw cseis_geolib::csException("Can't get first trace \n");

	/* Get parameters */
	if(!parObj.getparint("verbose",&verbose))	verbose = 1;
	if(!parObj.getparint("interp",&interp))	interp = 0;
	if(!parObj.getparfloat("value",&value))	value = 0.0;
        parObj.checkpars();

	/* Loop over traces */
	do{
		++itr;
      		for(i=0; i<tr.ns; ++i){
		    if(!isfinite(tr.data[i])) {
		       if (verbose)
	                warn("found NaN trace = %d  sample = %d", itr, i);

			if (interp) { /* interpolate nearest neighbors */
				      /* for NaN replacement value     */
				if (i==0 && isfinite(tr.data[i+1])) { 
					tr.data[i]=tr.data[i+1];
				} else if(i==tr.ns-1 && isfinite(tr.data[i-2])) {
					tr.data[i]= tr.data[i-2];
				} else if( isfinite(tr.data[i-1]) &&
						isfinite(tr.data[i+1]) ) {
					tr.data[i]=(tr.data[i-1]+tr.data[i+1])/2.0;
				}
			}
				
			/* use user defined NaNs replacement value */
            	       	tr.data[i] = value;
			}
		    }

      		su2cs->putTrace(&tr);
	} while(cs2su->getTrace(&tr));

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
