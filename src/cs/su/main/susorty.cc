/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSORTY: $Revision: 1.5 $ ; $Date: 2011/11/17 00:03:38 $	*/

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
std::string sdoc_susorty =
"								"
" SUSORTY - make a small 2-D common shot off-end  		"
"	    data set in which the data show geometry 		"
"	    values to help visualize data sorting.		"
"								"
"  susorty [optional parameters] > out_data_file  		"
"								"
" Optional parameters:						"
"	nt=100 		number of time samples			"
"	nshot=10 	number of shots				"
"	dshot=10 	shot interval (m)			"
"	noff=20 	number of offsets			"
"	doff=20 	offset increment (m)			"
"								"
" Notes:							"
" Creates a common shot su data file for sort visualization	"
"	       time samples           quantity			"
"	       ----------------      ----------			"
"	       first   25%           shot coord			"
"	       second  25%           rec coord			"
"	       third   25%           offset			"
"	       fourth  25%           cmp coord			"
"								"
"								"
" 1. default is shot ordered (hsv2 cmap looks best to me)	"
" susorty | suximage legend=1 units=meters cmap=hsv2		"
"								"
" 2. sort on cmp (note random order within a cmp)		"
" susorty | susort cdp > junk.su 				"
" suximage < junk.su legend=1 units=meters cmap=hsv2		"
"								"
" 3. sort to cmp and subsort on offset 	 			"
" susorty | susort cdp offset > junk.su 			"
" suximage < junk.su legend=1 units=meters cmap=hsv2		"
"								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace susorty {


/* Credits:
 *	CWP: Chris Liner  10.09.01
 *
 * Trace header fields set: ns, dt, sx, gx, offset, cdp, tracl 
 */
/**************** end self doc ***********************************/


segy tr;

void* main_susorty( void* args )
{
	int nt,nshot,noff;
        int ishot,ioff,it;
        float dt,dshot,doff,sx,gx,offset,cmp;

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_susorty );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	nt = 100;	parObj.getparint("nt", &nt);
	CHECK_NT("nt",nt);				tr.ns = nt;
	nshot = 10;	parObj.getparint("nshot", &nshot);
	noff  = 24;	parObj.getparint("noff", &noff);
	dt = 0.004;	parObj.getparfloat("dt", &dt);		tr.dt = dt*1000000;
	dshot = 10;	parObj.getparfloat("dshot", &dshot);
	doff = 20;	parObj.getparfloat("doff", &doff);
        parObj.checkpars();
	for (ishot = 0; ishot < nshot; ishot++) {
          sx = ishot*dshot;
	  for (ioff = 0; ioff < noff; ioff++) {
                offset = (ioff+1)*doff;
                gx = sx + offset; 
                cmp = (sx + gx)/2.;
		memset( (void *) tr.data, 0, nt * FSIZE);
                for (it = 0; it < nt/4; it++) {
		  tr.data[it] = sx;  
                }
                for (it = nt/4; it < nt/2; it++) {
		  tr.data[it] = gx;  
                }
                for (it = nt/2; it < 3*nt/4; it++) {
		  tr.data[it] = offset;  
                }
                for (it = 3*nt/4; it < nt; it++) {
		  tr.data[it] = cmp;  
                }
		tr.sx = sx;
		tr.gx = gx;
		tr.offset = offset;
		tr.cdp = cmp;
		tr.tracl = ishot*nshot + ioff + 1;
		su2cs->putTrace(&tr);
          }
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
