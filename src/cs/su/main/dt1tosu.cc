/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.		       */

/* DT1TOSU: $Revision: 1.17 $ ; $Date: 2012/03/30 20:35:28 $		*/

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
std::string sdoc_dt1tosu =
"								"
" DT1TOSU - Convert ground-penetrating radar data in the	"
"	Sensors & Software X.DT1 GPR format to SU format.	"
"								"
" dt1tosu < gpr_data_in_dt1_format  > stdout			"
"								"
" Optional parameters:						"
" ns=from header	number of samples per trace		"
" dt=.8		time sample interval (see below)		"
" swap=endian	endian is auto-determined =1 (big endian) swap	"
"		=0 don't swap bytes (little endian machines)	"
" verbose=0	silent						"
"		=1 S & S header values from first trace		"
"			sent to outpar				"
"		=2 S & S header values from all traces		"
"			sent to outpar				"
" outpar=/dev/tty	output parameter file			"
" list=0	silent						"
"		=1 list explaining labels used in verbose	"
"		     is printed to stderr			"
"								"
" Caution: An incorrect ns field will munge subsequent processing."
"								"
" Notes:							"
" For compatiblity with SEGY header, apparent dt is set to	"
" .8 ms (800 microsecs).  Actual dt is .8 nanosecs.		"
" Using TRUE DISTANCES, this scales velocity			"
" and frequency by a factor of 1 million.			"
"	Example: v_air = 9.83X10^8 ft/s	 (real)			"
"		 v_air = 983 ft/s	(apparent for su)	"
"	Example: fnyquist = 625 MHz	(real)			"
"		fnyquist = 625 Hz	(apparent for su)	"
"								"
" IBM RS6000, NeXT, SUN are examples of big endian machines	"
" PC's and DEC are examples of little endian machines		"
"								"
" Caveat:							"
" This program has not been tested on DEC, some modification of the"
" byte swapping routines may be required.			"
"								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace dt1tosu {


/* Credits:
 *	CWP: John Stockwell, Jan 1994   Based on a code "sugpr" by
 *	UTULSA: Chris Liner & Bill Underwood  (Dec93)
 * modifications permit S & S dt1 header information to be transferred
 * directly to SU header
 *
 * March 2012: CWP John Stockwell  updated for the revised
 * S&S DT1, which they still call "DT1" though it is different.
 *
 * Trace header fields set: ns, tracl, tracr, dt, delrt, trid,
 *			    hour, minute, second
 *
 * Reference: Sensors & Software pulseEKKO and Noggin^plus Data File
 *	     Formats
 * Publication of:
 * Sensors & Software: suburface imaging solutions
 * 1091 Brevik Place
 * Mississauga, ON L4W 3R7 Canada
 * Sensors & Software In
 * Tel: (905) 624-8909
 * Fax (905) 624-9365
 * E-mail: sales@sensoft.ca
 * Website: www.sensoft.ca
 */
/**************** end self doc *******************************************/

/* define hed structure */
#define SSHDRBYTES 128 /* size of ssdt1 */
typedef struct {
	float  tracl;	/* trace number */
	float  posit;	/* position */
	float  ns;	/* number of points per trace */
	float  topog;	/* Topographic data if any */
	float  nua;	/* not used "a" */ 
	float  bpp;	/* bytes/point (always 2 for Rev 3 firmware) */ 
	float  twind;	/* time window */ 
	float  nstks;	/* # of stacks */ 
	double  gps_x;	/* GPS X X-position */ 
	double  gps_y;	/* GPS Y Y-position */ 
	double  gps_z;	/* GPS Z Z-position */ 
	float  rx;	/* receiver x-position */ 
	float  ry;	/* receiver y-position */ 
	float  rz;	/* receiver z-position */ 
	float  tx;	/* transmitter x-position */ 
	float  ty;	/* transmitter y-position */ 
	float  tz;	/* transmitter z-position */ 
	float  tza;	/* time zero adjustment */ 
			/* where point(x) = point(x + adjustment) */
	float  zflag;	/* zero flag: 0=data ok ,  1 = zero data */
	float  nuc;	/* not used "c" */ 
	float  tod;	/* time of day data collected */
			/* in secs past midnight */ 
	float  cflag;	/* comment flag =1 if comment is attached */
	char comm[28];	/* comment 28 characters */
	short data[SU_NFLTS];  /* use SU maximum number data values */
} ssdt1;

ssdt1 sstr;
segy tr;

/* list explaining the ssdt1 convention */
char *list[] = {
" float	 tracl;		trace number				",
" float	 posit;		position				",
" float	 ns;		number of points per trace		",
" float	 topog;		Topographic data if any			",
" float	 nua;		not used \"a\"				",
" float	 bpp;		bytes/point (always 2 for Rev 3 firmware)",
" float	 twind;		time window				",
" float	 nstks;		number of stacks			",
" double  gps_x;	/* GPS X X-position */			",
" double  gps_y;	/* GPS Y Y-position */ 			",
" double  gps_z;	/* GPS Z Z-position */			", 
" float  rx;		/* receiver x-position */		", 
" float  ry;		/* receiver y-position */		",
" float  rz;		/* receiver z-position */		",
" float  tx;		/* transmitter x-position */		",
" float  ty;		/* transmitter y-position */ 		",
" float  tz;		/* transmitter z-position */ 		",
" float	 tza;		time zero adjustment			",
"			where point(x) = point(x + adjustment)	",
" float	 zflag;		zero flag: 0=data ok ,	1 = zero data	",
" float	 nuc;		not used \"c\"				",
" float	 tod;		time of day data are collected		",
"				in secs past midnight		",
" float	 cflag;		comment flag =1 if comment is attached	",
" char comm[28];	comment 28 characters			",
" short data[SU_NFLTS];	data part of the ssdt1 structure	",
NULL};

/* pointer to list */
char **listptr=list;

/* function prototypes for internally defined subroutines */
void swap_ss_sstr(void);
void fprintf_sstr(FILE *outparfp);


void* main_dt1tosu( void* args )
{
	int i,ns;		/* counter, number of samples */
	float dt;		/* time sample interval */	
	size_t ssdatabytes;	/* number of bytes in SS data portion */
	size_t ssbytes;		/* total number of SS bytes */
	int swap,verbose,list;	/* flags */
	float hour,minute;	/* hour and minute */
	char *outpar=NULL;	/* outpar filename */
	FILE *outparfp=NULL;	/* outpar file pointer */

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_dt1tosu );
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
	default: /* rest are OK */
	break;
	}

	/* Determin endian of system */
	if (!parObj.getparint("swap", &swap))      {
	      union { short s; char c[2]; } testend;
		testend.s = 1;
 		swap = (testend.c[0] == '\0') ? 1 : 0;
	}

	/* Read header information from first header */
	efread(&sstr, 1, SSHDRBYTES, stdin);
	if (swap) swap_ss_sstr();

	/* Get parameters */
	if (!parObj.getparint("ns", &ns)) 	 ns = sstr.ns ;
	if (!ns) throw cseis_geolib::csException("ns not set!");
	warn("ns= %d", ns);

	if (!parObj.getparfloat("dt", &dt))		dt= .8 ;
	dt *= 1000.;
	if (!parObj.getparint("list", &list))		list= 0;
	if (!parObj.getparint("verbose", &verbose))	verbose= 0;
	if (!parObj.getparstring("outpar", &outpar))	outpar = "/dev/tty" ;
	parObj.checkpars();

	outparfp = efopen(outpar, "w");

	/* print out ss header field for first trace */
	if (verbose==1) fprintf_sstr(outparfp);
	if (list) while (*listptr) fprintf(stderr,"%s\n", *listptr++);

	/* assign tr header values */
	tr.tracl = sstr.tracl;

	tr.sx = (int) sstr.tx;
	tr.sy = (int) sstr.ty;
	tr.selev = (int) sstr.tz;

	tr.gx = (int) sstr.rx;
	tr.gy = (int) sstr.ry;
	tr.gelev = (int) sstr.rz;

	tr.ns = sstr.ns;	
	tr.dt = NINT(dt); 
	tr.delrt = sstr.tza;
	tr.trid = sstr.zflag + 1;

	/* compute hour and minute on first trace */
	hour = sstr.tod/3600;
	minute = (hour - (int) hour) * 60;

	/* set hour, minute, and second on first trace */
	tr.hour = hour;
	tr.minute = minute;
	tr.sec = ( minute - (int) minute ) * 60;

	
	/* read in and swap the data portion of the first trace */
	/* loop over data values in X.dt1 data */
	for (i=0 ; i < ns ; i++) {
		short temp;	/* temporary variable */
				/* to store data vals in */
		efread(&temp, 1, sizeof(short), stdin);

		if (swap) swap_short_2(&temp); 

		tr.data[i] = (float) temp;
	}

	su2cs->putTrace(&tr);

	/* assign size of data portion on first trace */
	ssdatabytes = ns*sizeof(short);
	ssbytes = SSHDRBYTES + ssdatabytes;

	/* load selected sstr values and data into tr */
	while (efread(&sstr,1,ssbytes,stdin)) {

		/* swap header values */
		if (swap) swap_ss_sstr();

		/* assign tr header values */
		tr.tracl = sstr.tracl;

		tr.sx = (int) sstr.tx;
		tr.sy = (int) sstr.ty;
		tr.selev = (int) sstr.tz;

		tr.gx = (int) sstr.rx;
		tr.gy = (int) sstr.ry;
		tr.gelev = (int) sstr.rz;

		tr.ns = sstr.ns;	
		tr.dt = NINT(dt); 
		tr.delrt = sstr.tza;
		tr.trid = sstr.zflag + 1;

		/* compute hour and minute on first trace */
		hour = sstr.tod/3600;
		minute = (hour - (int) hour) * 60;

		/* set hour, minute, and second on first trace */
		tr.hour = hour;
		tr.minute = minute;
		tr.sec = ( minute - (int) minute ) * 60;

		/* loop over data values in X.dt1 data */
		for (i=0 ; i < ns ; i++) {
			short temp;	/* temporary variable */
				/* to store data vals in */

			temp = sstr.data[i];
			if (swap) swap_short_2(&temp); 

			tr.data[i] = (float) temp;
		}
		if (verbose==2) fprintf_sstr(outparfp);
	
		su2cs->putTrace(&tr);
	}


	if (verbose) efclose(outparfp);
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


/* swap */

void swap_ss_sstr(void)
{
	swap_float_4(&sstr.tracl);
	swap_float_4(&sstr.posit);
	swap_float_4(&sstr.ns);
	swap_float_4(&sstr.topog);
	swap_float_4(&sstr.bpp);
	swap_float_4(&sstr.twind);
	swap_float_4(&sstr.nstks);
	swap_double_8(&sstr.gps_x);
	swap_double_8(&sstr.gps_y);
	swap_double_8(&sstr.gps_z);
	swap_float_4(&sstr.rx);
	swap_float_4(&sstr.ry);
	swap_float_4(&sstr.rz);
	swap_float_4(&sstr.tx);
	swap_float_4(&sstr.ty);
	swap_float_4(&sstr.tz);
	swap_float_4(&sstr.tza);
	swap_float_4(&sstr.zflag);
	swap_float_4(&sstr.nuc);
	swap_float_4(&sstr.tod);
}



void fprintf_sstr(FILE *outparfp)
{ /* send ssdt1 values to outpar */
	int tracl=sstr.tracl; /* first trace counter */
	float  posit=sstr.posit; /* position */
	int ns=sstr.ns;	 /* number of samples */
	float  topog=sstr.topog; /* topographic data */
	int bpp=sstr.bpp;	/* bytes per point */ 
	float twind=sstr.twind;	/* time window */
	int nstks=sstr.nstks;	/* number of stacks */
	float  tza=sstr.tza;	/* time zero adjustment */
	int  zflag=sstr.zflag;	/* zero flag */
	float nuc=sstr.nuc; /* unused "c" */
	int  cflag=sstr.cflag; /* comment flag */


	fprintf(outparfp, "tracl = %d\n", tracl);
	fprintf(outparfp, "posit = %f\n", posit);
	fprintf(outparfp, "ns = %d\n", ns);
	fprintf(outparfp, "topog = %f\n", topog);
	fprintf(outparfp, "bpp = %d\n", bpp);
	fprintf(outparfp, "twind = %f\n", twind);
	fprintf(outparfp, "nstks = %d\n", nstks);
	fprintf(outparfp, "tza = %f\n", tza);
	fprintf(outparfp, "zflag = %d\n", zflag);
	fprintf(outparfp, "nuc = %f\n", nuc);
	fprintf(outparfp, "cflag = %d\n", cflag);
	if(cflag) efwrite(sstr.comm,28,sizeof(char),outparfp);
	else      fprintf(outparfp, "no comm\n");
}

} // END namespace
