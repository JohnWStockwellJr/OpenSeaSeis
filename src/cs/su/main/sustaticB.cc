/* Copyright (c) Colorado School of Mines, 2010.*/
/* All rights reserved.                       */

/* SUSTATIC: $Revision: 1.1 $ ; $Date: 2012/07/19 18:13:34 $	*/

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

/*********************** self documentation ******************************/
std::string sdoc_sustaticB =
"									"
" SUSTATICB - Elevation static corrections, apply corrections from	"
"	      headers or from a source and receiver statics file	"
"	      (beta submitted by J. W. Neese)				"
"									"
"     sustaticB <stdin >stdout  [optional parameters]	 		"
"									"
" Required parameters:							"
"	none								"
" Optional Parameters:							"
"	v0=v1 or user-defined	or from header, weathering velocity	"
"	v1=user-defined		or from header, subweathering velocity	"
"	hdrs=0			=1 to read statics from headers		"
" 				=2 to read statics from files		"
"				=3 to read from output files of suresstat"
"	sign=1			apply static correction (add tstat values)"
"				=-1 apply negative of tstat values	"
" Options when hdrs=2 and hdrs=3:					"
"	sou_file=		input file for source statics (ms) 	"
"	rec_file=		input file for receiver statics (ms) 	"
"	ns=240 		(2)number of sources; (3) max fldr	"
"	nr=335 			number of receivers 			"
"	no=96 			number of offsets			"
"									"
" Notes:								"
" For hdrs=1, statics calculation is not performed, statics correction  "
" is applied to the data by reading statics (in ms) from the header.	"
"									"
" For hdrs=0, field statics are calculated, and				"
" 	input field sut is assumed measured in ms.			"
" 	output field sstat = 10^scalel*(sdel - selev + sdepth)/swevel	"
" 	output field gstat = sstat - sut/1000.				"
" 	output field tstat = sstat + gstat + 10^scalel*(selev - gelev)/wevel"
"									"
" For hdrs=2, statics are surface consistently obtained from the 	"
" statics files. The geometry should be regular.			"
" The source- and receiver-statics files should be unformated C binary 	"
" floats and contain the statics (in ms) as a function of surface location."
"									"
" For hdrs=3, statics are read from the output files of suresstat, with "
" the same options as hdrs=2 (but use no=max traces per shot and assume "
" that ns=max fldr number and nr=max receiver number).			"
" For each shot number (trace header fldr) and each receiver number     "
" (trace header tracf) the program will look up the appropriate static  "
" correction.  The geometry need not be regular as each trace is treated"
" independently.							"
"									"
" Caveat:  The static shifts are computed with the assumption that the  "
" desired datum is sea level (elevation=0). You may need to shift the	"
" selev and gelev header values via  suchw.				"
" Example: subtracting min(selev,gelev)=25094431			"
"									"
" suchw < CR290.su key1=selev,gelev key2=selev,gelev key3=selev,gelev \\ "
"            a=-25094431,-25094431 b=1,1 c=0,0 > CR290datum.su		"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sustaticB {


/* Credits:
 *	CWP: Jamie Burns
 *
 *	CWP: Modified by Mohammed Alfaraj, 11/10/1992, for reading
 *	     statics from headers and including sign (+-) option
 *
 *      CWP: Modified by Timo Tjan, 29 June 1995, to include input of
 *           source and receiver statics from files. 
 *
 *	modified by Thomas Pratt, USGS, Feb, 2000 to read statics from
 * 	     the output files of suresstat
 *
 * Logic changed by JWN to fix options hdrs=2,3 ???
 *
 * Trace header fields accessed:  ns, dt, delrt, gelev, selev,
 *	sdepth, gdel, sdel, swevel, sut, scalel, fldr, tracf
 * Trace header fields modified:  sstat, gstat, tstat
 */

/************************ end self doc ***********************************/


segy intrace, outtrace;

void* main_sustaticB( void* args )
{
	int nt;		/* number of samples on output trace	*/
	float dt;	/* sample rate on outpu trace		*/
	int itime;	/* counter          			*/
	float tmin;	/* first time sample on output trace	*/
	float tsd=0.0;	/* time to move source to datum         */
	float trd=0.0;	/* time to move 0 offset receiver       */
	float v0;	/* weathering velocity			*/
	float v1;	/* subweathering velocity		*/
	int hdrs; 	/* flag to read statics from headers	*/ 
	float *t;	/* array of output times		*/
	float tstat;	/* total (source and receiver) statics	*/
	int sign;	/* to add (+) or subtract (-) statics	*/
	int no;		/* number of offsets per shot 		*/
	int io;		/* offset counter 			*/
	int is;		/* source counter 			*/
	int ir;		/* receiver counter 			*/
	int ns;		/* number of sources = number of source statics */
	int nr;		/* number of receiver = number of rec. statics	*/
	float *sou_statics=NULL;	/* array of source statics	*/
	float *rec_statics=NULL;	/* array of receiver statics	*/
	FILE *fps, *fpr;	/* file pointers for statics input 	*/
	cwp_String sou_file, rec_file; /* statics filenames 		*/

	/* Hook up getpar */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sustaticB );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get information from first trace */
	if (!cs2su->getTrace(&intrace)) throw cseis_geolib::csException("can't get first trace");
	nt   = intrace.ns;
	tmin = intrace.delrt/1000.0;
	dt   = ((double) intrace.dt)/1000000.0;
	
	/* Get parameters */
	if (!parObj.getparfloat("v1", &v1))          v1 = (float) intrace.swevel;
	if (!parObj.getparfloat("v0", &v0))
                v0 = (float) ((intrace.wevel) ? intrace.wevel : v1);
	if (!parObj.getparint("hdrs", &hdrs))        hdrs = 0;
	if (!parObj.getparint("sign", &sign))        sign = 1;

	/* Allocate vector of output times */
	t = ealloc1float(nt);

	/* reading source and receiver statics from files */
	if ((hdrs == 2) || (hdrs == 3)){

		/* parObj.getpar statics file related parameters */
		if (!parObj.getparint("ns", &ns))        ns = 240;
		if (!parObj.getparint("nr", &nr))        nr = 335;
		if (!parObj.getparint("no", &no))        no = 96;

		/* parObj.getpar statics file names */
        	parObj.getparstring("sou_file",&sou_file);
        	parObj.getparstring("rec_file",&rec_file);

		/* allocate space */
		rec_statics = alloc1float(nr);
        	sou_statics = alloc1float(ns);

		/* open and read from receiver statics file */
        	if((fpr=efopen(rec_file,"rb"))==NULL)
                	throw cseis_geolib::csException("cannot open stat_file=%s\n",rec_file);
        	efread(rec_statics, sizeof(float),nr,fpr);
        	efclose(fpr);

		/* open and read from source statics file */
        	if((fps=efopen(sou_file,"rb"))==NULL)
                	throw cseis_geolib::csException("cannot open stat_file=%s\n",sou_file);
        	efread(sou_statics, sizeof(float),ns,fps);
        	efclose(fps);
	}

	/* Initialize tstat */
	tstat = 0.0;

	/* Loop on traces */	
	io = 0; is = 0;
	do {
		int temp = SGN(intrace.scalel)*log10(abs((int) intrace.scalel));
		float scale;
                scale = pow(10., (float)temp);
		
		/* copy and adjust header */
		memcpy( (void *) &outtrace, (const void *) &intrace, HDRBYTES);
	
		/* compute static correction if necessary */
		if(!hdrs) {
		    	tsd = scale *
			(-intrace.selev + intrace.sdel + intrace.sdepth)/v1;
			trd = tsd - intrace.sut/1000.0;
			tstat = tsd + trd +
				scale * (intrace.selev - intrace.gelev)/v0;

		/* else, read statics from headers */
		} else { 
			/* Initialize header field for output trace */
			outtrace.sstat = intrace.sstat;
			outtrace.gstat = intrace.gstat;
			outtrace.tstat = intrace.sstat+intrace.gstat;
			if (hdrs == 1) {
				tstat = outtrace.tstat/1000.0;
			}
			if (hdrs == 2) {
				ir = is + io;
				if (is <= ns) tsd = sou_statics[is]/1000.0;
				if (ir > 0 && ir <= nr)
					trd = rec_statics[ir]/1000.0;

				tstat = tsd + trd;
				io ++;
				if (io > no-1) {
					io = 0; is++;
				}
			}
			if (hdrs == 3) {
				tsd = sou_statics[intrace.fldr]/1000.0;
				trd = rec_statics[intrace.tracf]/1000.0;

				tstat = tsd + trd;
			}
		}
		
		/* Compute output times */
		for (itime=0; itime<nt; ++itime)
			t[itime] = tmin + itime*dt + sign*tstat;

		/* sinc interpolate new data */
		ints8r(nt, dt, tmin, intrace.data, 
				0.0, 0.0, nt, t, outtrace.data);
		
		/* set header field for output trace */
		if(hdrs == 0 || hdrs == 2 || hdrs == 3){

			/* value is added to existing header values */
			/* this permits multiple static corrections */
			outtrace.sstat += (1000.0 * tsd);
			outtrace.gstat += (1000.0 * trd);
			outtrace.tstat += (1000.0 * tstat);
		} 
		
		su2cs->putTrace(&outtrace);
	} while (cs2su->getTrace(&intrace));


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
