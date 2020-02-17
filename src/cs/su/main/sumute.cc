/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUMUTE: $Revision: 1.39 $ ; $Date: 2011/11/17 00:03:38 $	*/

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
std::string sdoc_sumute =
" 	   								"
" SUMUTE - MUTE above (or below) a user-defined polygonal curve with	" 
"	   the distance along the curve specified by key header word 	"
" 	   								"
" sumute <stdin >stdout xmute= tmute= [optional parameters]		"
" 									"
" Required parameters:							"
" xmute=		array of position values as specified by	"
" 			the `key' parameter				"
" tmute=		array of corresponding time values (sec)	"
" 			in case of air wave muting, correspond to 	"
" 			air blast duration				"
"  ... or input via files:						"
" nmute=		number of x,t values defining mute		"
" xfile=		file containing position values as specified by	"
" 			the `key' parameter				"
" tfile=		file containing corresponding time values (sec)	"
"  ... or via header:							"
" hmute=		key header word specifying mute time		"
" 									"
" Optional parameters:							"
" key=offset		Key header word specifying trace offset 	"
" 				=tracl  use trace number instead	"
" ntaper=0		number of points to taper before hard		"
"			mute (sine squared taper)			"
" mode=0	   mute ABOVE the polygonal curve			"
"		=1 to zero BELOW the polygonal curve			"
"		=2 to mute below AND above a straight line. In this case"
"		 	xmute,tmute describe the total time length of   "
"			the muted zone as a function of xmute the slope "
"			of the line is given by the velocity linvel=	"
"	 	=3 to mute below AND above a constant velocity hyperbola"
"			as in mode=2 xmute,tmute describe the total time"
"			length of the mute zone as a function of xmute,  "
"			the velocity is given by the value of linvel=	"
" 		=4 to mute below AND above a user defined polygonal line"
"			given by xmute, tmute pairs. The widths in time " 
"			of the muted zone are given by the twindow vector"
" linvel=330   		constant velocity for linear or hyperbolic mute	"
" tm0=0   		time shift of linear or hyperbolic mute at	"
"			 \'key\'=0					"
" twindow=	vector of mute zone widths in time, operative only in mode=4"
"  ... or input via file:						"
" twfile= 								"
"									"
" Notes: 								"
" The tmute interpolant is extrapolated to the left by the smallest time"
" sample on the trace and to the right by the last value given in the	"
" tmute array.								"
"									"
" The files tfile and xfile are files of binary (C-style) floats.	"
"									"
" In the context of this program \"above\" means earlier time and	"
" \"below\" means later time (above and below as seen on a seismic section."
"									"
" The mode=2 option is intended for removing air waves. The mute is	"
" is over a narrow window above and below the line specified by the	"
" the line specified by the velocity \"linvel\". Here the values of     "
" tmute, xmute or tfile and xfile define the total time width of the mute."
"									"
" If data are spatial, such as the (z-x) output of a migration, then    "
" depth values are used in place of times in tmute and tfile. The value "
" of the depth sampling interval is given by the d1 header field	"
"									"
" Caveat: if data are seismic time sections, then tr.dt must be set. If "
" data are seismic depth sections, then tr.trid must be set to the value"
" of TRID_DEPTH and the tr.d1 header field must be set.			"
" To find the value of TRID_DEPTH:  					"
" type: 								"
"     sukeyword trid							"
"	and look for the entry for \"Depth-Range (z-x) traces\"		"
"  									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sumute {


/* Credits:
 *
 *	SEP: Shuki Ronen
 *	CWP: Jack K. Cohen, Dave Hale, John Stockwell
 *	DELPHI: Alexander Koek added airwave mute.
 *      CWP: John Stockwell added modes 3 and 4
 *	USBG: Florian Bleibinhaus added hmute + some range checks on mute times
 * Trace header fields accessed: ns, dt, delrt, key=keyword, trid, d1
 * Trace header fields modified: muts or mute
 */
/**************** end self doc ***********************************/

segy tr;


#define SQ(x) ((x))*((x))

void* main_sumute( void* args )
{
	char *key=NULL;		/* header key word from segy.h		*/	
	char *type=NULL;	/* ... its type				*/
	int index;		/* ... its index			*/
	Value val;		/* ... its value			*/
	float fval;		/* ... its value cast to float		*/
	float twfval;		/* ... its value cast to float		*/
	float *xmute=NULL;	/* array of key mute curve values	*/
	float *tmute=NULL;	/* ...		mute curve time values 	*/
	float *twindow=NULL;	/* ...	mute window time values mode=4	*/
	float linvel;		/* linear velocity			*/
	float tm0;		/* time shift of mute=2 or 3 for 'key'=0*/
	float *taper=NULL;	/* ...		taper values		*/
	int nxmute=0;		/* number of key mute values		*/
	int ntmute;	/* ...		mute time values 	*/
	int ntwindow;	/* ...		mute time values 	*/
	int ntaper;	/* ...		taper values		*/

	int below;	/* mute below curve			*/
	int mode;	/* kind of mute (top, bottom, linear)	*/
	int absolute;    /* Take absolute value of key for mode=2 */
	int hmute=0;	/* read mute times from header		*/

	int nxtmute;	/* number of mute values 		*/
	cwp_String xfile="";	/* file containing positions by key	*/
	FILE *xfilep;		/* ... its file pointer			*/
	cwp_String tfile="";	/* file containing times	 	*/
	FILE *tfilep;		/* ... its file pointer			*/

	cwp_String twfile="";	/* file containing mute time windows 	*/
	FILE *twfilep;		/* ... its file pointer			*/

	cwp_Bool seismic;	/* cwp_true if seismic, cwp_false not seismic */

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sumute );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get parameters */
	if (!parObj.getparint("mode", &mode))		mode = 0;

	if (parObj.getparstring("hmute", &key)) { hmute = 1; } else

	if (!(parObj.getparstring("xfile",&xfile) && getparstring("tfile",&tfile))) {
		if (!(nxmute = parObj.countparval("xmute")))
			throw cseis_geolib::csException("must give xmute= vector");
		if (!(ntmute = parObj.countparval("tmute")))
			throw cseis_geolib::csException("must give tmute= vector");
		if (nxmute != ntmute)
			throw cseis_geolib::csException("lengths of xmute, tmute must be the same");
		xmute = ealloc1float(nxmute);	parObj.getparfloat("xmute", xmute);
		tmute = ealloc1float(nxmute);	parObj.getparfloat("tmute", tmute);
		if (mode==4) {
			if (!(ntwindow = parObj.countparval("twindow")))
				throw cseis_geolib::csException("must give twindow= vector");
			if (nxmute != ntwindow)
				throw cseis_geolib::csException("lengths of xmute, twindow must be the same");
			twindow = ealloc1float(nxmute);
			parObj.getparfloat("twindow", twindow);
		}
	} else {
		CSMUSTGETPARINT("nmute",&nxtmute);
		nxmute = nxtmute;
		xmute = ealloc1float(nxtmute);
		tmute = ealloc1float(nxtmute);

		if((xfilep=fopen(xfile,"r"))==NULL)
			throw cseis_geolib::csException("cannot open xfile=%s\n",xfile);
		if (fread(xmute,sizeof(float),nxtmute,xfilep)!=nxtmute)
			throw cseis_geolib::csException("error reading xfile=%s\n",xfile);
		fclose(xfilep);

		if((tfilep=fopen(tfile,"r"))==NULL)
			throw cseis_geolib::csException("cannot open tfile=%s\n",tfile);
		if (fread(tmute,sizeof(float),nxtmute,tfilep)!=nxtmute)
			throw cseis_geolib::csException("error reading tfile=%s\n",tfile);
		fclose(tfilep);

		if (mode==4) {
			if((twfilep=fopen(twfile,"r"))==NULL)
				throw cseis_geolib::csException("cannot open tfile=%s\n",twfile);
			if (fread(twindow,sizeof(float),nxtmute,twfilep)!=nxtmute)
				throw cseis_geolib::csException("error reading tfile=%s\n",tfile);
			fclose(twfilep);
		}
	}

	if (!parObj.getparint("ntaper", &ntaper))	ntaper = 0;
	if (parObj.getparint("below", &below))	{
		mode = below; 
		warn ("use of below parameter is obsolete. mode value set to %d \n", mode);
	}
	if (!parObj.getparint("absolute", &absolute))		absolute = 1;
	if (hmute==0) if (!parObj.getparstring("key", &key))	key = "offset";
	if (!parObj.getparfloat("linvel", &linvel))	linvel = 330;
	if (!parObj.getparfloat("tm0", &tm0))		tm0 = 0;
        parObj.checkpars();
	if (linvel==0) throw cseis_geolib::csException("linear velocity can't be 0");

	/* get key type and index */
	type = hdtype(key);
	index = getindex(key);

	/* Set up taper weights if tapering requested */
	if (ntaper) {
		register int k;
		taper = ealloc1float(ntaper);
		for (k = 0; k < ntaper; ++k) {
			float s = sin((k+1)*PI/(2*ntaper));
			taper[k] = s*s;
		}
	}

						
	/* Get info from first trace */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't read first trace");

	seismic = ISSEISMIC(tr.trid);

	if (seismic) {
		if (!tr.dt) throw cseis_geolib::csException("dt header field must be set");
	} else if (tr.trid==TRID_DEPTH) {   /* depth section */
		if (!tr.d1) throw cseis_geolib::csException("d1 header field must be set");
	} else {
		throw cseis_geolib::csException("tr.trid = %d, unsupported trace id",tr.trid);
	}

	/* Loop over traces */
	do {
		int nt     = (int) tr.ns;
		float tmin = tr.delrt/1000.0;
		float dt = ((double) tr.dt)/1000000.0;
		float t;
		float tw;
		int nmute;
		int itaper;
		int topmute;
		int botmute;
		int ntair=0;
		register int i;

		if (!seismic) { 
			tmin = 0.0;
			dt = tr.d1;	
		}

		/* get value of key and convert to float */
		gethval(&tr, index, &val);
		fval = vtof(type,val);
		

		if (hmute==1) {
			t = fval/1000.;
		} else {
		/* linearly interpolate between (xmute,tmute) values */
		intlin(nxmute,xmute,tmute,tmin,tmute[nxmute-1],1,&fval,&t); 
		}

		if (absolute) fval = abs(fval);
		/* do the mute */
		if (mode==0) {	/* mute above */
			nmute = MIN(NINT((t - tmin)/dt),nt);
			if (nmute>0) memset( (void *) tr.data, 0, nmute*FSIZE);
			for (i = 0; i < ntaper; ++i)
				if (i+nmute>0) tr.data[i+nmute] *= taper[i];
			if (seismic) {
				tr.muts = NINT(t*1000);
			} else	{
				tr.muts = NINT(t);
			}
		} else if (mode==1){	/* mute below */
			nmute = MAX(0,NINT((tmin + nt*dt - t)/dt));
			memset( (void *) (tr.data+nt-nmute), 0, nmute*FSIZE);
			for (i = 0; i < ntaper; ++i)
				if (nt>nmute+i && nmute+i>0)
					tr.data[nt-nmute-1-i] *= taper[i];
			if (seismic) {
				tr.mute = NINT(t*1000);
			} else	{
				tr.mute = NINT(t);
			}
		} else if (mode==2){	/* air wave mute */
			nmute = NINT((tmin+t)/dt);
			ntair=NINT(tm0/dt+fval/linvel/dt);
			topmute=MIN(MAX(0,ntair-nmute/2),nt);
			botmute=MIN(nt,ntair+nmute/2);
			memset( (void *) (tr.data+topmute), 0,
					(botmute-topmute)*FSIZE);
			for (i = 0; i < ntaper; ++i){
				itaper=ntair-nmute/2-i;
				if (itaper > 0) tr.data[itaper] *=taper[i];
			}	
			for (i = 0; i < ntaper; ++i){
				itaper=ntair+nmute/2+i;
				if (itaper<nt) tr.data[itaper] *=taper[i];
			}
		} else if (mode==3) {	/* hyperbolic mute */
			nmute = NINT((tmin + t)/dt);
			ntair=NINT(sqrt( SQ((float)(tm0/dt))+SQ((float)(fval/linvel/dt)) ));
			topmute=MIN(MAX(0,ntair-nmute/2),nt);
			botmute=MIN(nt,ntair+nmute/2);
			memset( (void *) (tr.data+topmute), 0,
					(botmute-topmute)*FSIZE);
			for (i = 0; i < ntaper; ++i){
				itaper=ntair-nmute/2-i;
				if (itaper > 0) tr.data[itaper] *=taper[i];
			}	
			for (i = 0; i < ntaper; ++i){
				itaper=ntair+nmute/2+i;
				if (itaper<nt) tr.data[itaper] *=taper[i];
			}
		} else if (mode==4) {	/* polygonal mute */
			tmin=twindow[0];
			intlin(nxmute,xmute,twindow,tmin,twindow[nxmute-1],1,&twfval,&tw); 
			if (absolute) twfval = abs(twfval);

			nmute = NINT(tw/dt);
			ntair = NINT(t/dt);

			topmute=MIN(MAX(0,ntair-nmute/2),nt);
			botmute=MIN(nt,ntair+nmute/2);
			memset( (void *) (tr.data+topmute), 0,
					(botmute-topmute)*FSIZE);
			for (i = 0;i < ntaper; ++i){
				itaper=ntair-nmute/2-i;
				if (itaper > 0) tr.data[itaper] *=taper[i];
			}	
			for (i = 0; i < ntaper; ++i){
				itaper=ntair+nmute/2+i;
				if (itaper<nt) tr.data[itaper] *=taper[i];
			}


		}
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
