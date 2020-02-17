/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUVIBRO: $Revision: 1.14 $ ; $Date: 2012/10/22 14:35:30 $	*/

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

#define TWOPI 2.0*PI

/*********************** self documentation **********************/
std::string sdoc_suvibro =
"								"
" SUVIBRO - Generates a Vibroseis sweep (linear, linear-segment,"
"			dB per Octave, dB per Hertz, T-power)	"
"								"
" suvibro [optional parameters] > out_data_file			"
"								"
" Optional Parameters:						"
" dt=0.004		time sampling interval			"
" sweep=1	  	linear sweep			  	"
" 		  	=2 linear-segment			"
" 		  	=3 decibel per octave	 		"
" 		  	=4 decibel per hertz	  		"
" 		  	=5 t-power				"
" swconst=0.0		sweep constant (see note)		"
" f1=10.0		sweep frequency at start		"
" f2=60.0		sweep frequency at end			"
" tv=10.0		sweep length				"
" phz=0.0		initial phase (radians=1 default)	"
" radians=1		=0 degrees				"
" fseg=10.0,60.0	frequency segments (see notes)		"
" tseg=0.0,10.0		time segments (see notes)		"
" t1=1.0		length of taper at start (see notes)	"
" t2=1.0		length of taper at end (see notes)	"
" taper=1		linear					"
"		  	=2 sine					"
"			=3 cosine				"
"			=4 gaussian (+/-3.8)			"
"			=5 gaussian (+/-2.0)			"
"								"
" Notes:							"
" The default tapers are linear envelopes. To eliminate the	"
" taper, choose t1=t2=0.0.					"
"								"
" \"swconst\" is active only with nonlinear sweeps, i.e. when	"
" sweep=3,4,5.							" 
" \"tseg\" and \"fseg\" arrays are used when only sweep=2	"
"								"
" Sweep is a modulated cosine function.				"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suvibro {


/*
 * Author: CWP: Michel Dietrich
   Rewrite: Tagir Galikeev, CWP,  7 October 1994
 *
 * Trace header fields set: ns, dt, tracl, sfs, sfe, slen, styp
 */
/**************** end self doc ***********************************/

segy tr;

/* Prototypes for functions used interally */
void Linear( float fs, float fe, float T, float dt, float phz );
void Linear_Segment( float *freq, float *time, int isegm, float T, 
			float dt, float phz );
void dB_per_octave( float fs, float fe, float T, float dt, float swconst, 
			float phz );
void dB_per_hertz( float fs, float fe, float T, float dt, float swconst,
			float phz );
void t_power( float fs, float fe, float T, float dt, float swconst,
			float phz );
void sweep_taper( float t1, float t2, int tap_type, float T, float dt );


void* main_suvibro( void* args )
{
	/* sweep frequency parameters */
	float f1;		/* sweep frequency at beginning */
	float f2;		/* sweep frequency at end */
	float tv;		/* total sweep band */

	/* initial sweep phase paramters */
	int radians=1;		/* =1 phase is in radians, =0 degrees */
	float phz;		/* initial phase */

	/* tapering parameters */
	float t1; 		/* taper at beginning of sweep */
	float t2;		/* taper at end of sweep */

	float dt;		/* time sampling interval on output */
	float swconst;		/* sweep constant for nonlinear sweeps */

	int sweep, i;
	int nt, isegm=0, ntseg, nfseg, taper;
	float *fseg=NULL, *tseg=NULL, *time=NULL;


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suvibro );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	
	/* set parameters and fill header fields */
	if (!parObj.getparint("sweep", &sweep)) sweep = 1;
	if (!parObj.getparfloat("swconst", &swconst))	swconst = 0.0;
	if (!parObj.getparfloat("f1", &f1))	f1 = 10.0;
	if (!parObj.getparfloat("f2", &f2))	f2 = 60.0;
	if (!parObj.getparfloat("tv", &tv))	tv = 10.0;
	if (!parObj.getparint("radians", &radians)) radians = 1;
	if (!parObj.getparfloat("phz", &phz))	phz = 0.0;
	if (!parObj.getparfloat("t1", &t1))	t1 = 1.0;
	if (!parObj.getparfloat("t2", &t2))	t2 = 1.0;
	if (!parObj.getparint("taper", &taper)) taper = 1;
	if (!parObj.getparfloat("dt", &dt))	dt = 0.004;
	
	tr.dt = (float) (1000000.0*dt);
	if ((nt = tv/dt + 1) >= SU_NFLTS) throw cseis_geolib::csException("nt=tv/dt=%d -- too big", nt);
	if (t1 + t2 > tv)
		throw cseis_geolib::csException("sum of tapers t1=%f, t2=%f exceeds tv=%f", t1,t2,tv);

	/* initial phase */
	if(!radians) { /* then degrees */
		phz = phz*TWOPI/360.0;
	}
		
	if (sweep==2) {
		ntseg = parObj.countparval("tseg");
		if (ntseg==0) ntseg = 2;
		tseg = ealloc1float(ntseg);
		if (!parObj.getparfloat("tseg",tseg)) { tseg[0] = 0.0; tseg[1] = tv; }
		nfseg = parObj.countparval("fseg");
		if (nfseg==0) nfseg = 2;
		if (nfseg!=ntseg) throw cseis_geolib::csException("number of tseg and fseg must be equal");
		fseg = ealloc1float(nfseg);
		if (!parObj.getparfloat("fseg",fseg)) { fseg[0] = f1; fseg[1] = f2; }
		if (ntseg > 2)
		for (i=1; i<ntseg; ++i)
		if (tseg[i]<=tseg[i-1])
		  throw cseis_geolib::csException("tseg must increase monotonically");

		/* prepare freq and time arrays */
		isegm = ntseg - 1;
		tv = tseg[isegm];
		time = ealloc1float(ntseg-1);

		/* compute each segment length */
		for (i=0; i<isegm; i++)
		time[i] = tseg[i+1] - tseg[i];
	  }

	tr.ns = nt;
	tr.tracl = 1;
	tr.sfs = f1;
	tr.sfe = f2;
	tr.slen = 1000.0 * tv;
	tr.styp = sweep;		/* sweep id code */

	if (sweep==1) Linear( f1, f2, tv, dt, phz ); 
	if (sweep==2) Linear_Segment( fseg, time, isegm, tv, dt, phz );
	if (sweep==3) dB_per_octave( f1, f2, tv, dt, swconst, phz );
	if (sweep==4) dB_per_hertz( f1, f2, tv, dt, swconst, phz );
	if (sweep==5) dB_per_hertz( f1, f2, tv, dt, swconst, phz );

	if ( (t1!=0.) || (t2!=0.) )
	  sweep_taper( t1, t2, taper, tv, dt );

	su2cs->putTrace(&tr); 	 	
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


void Linear ( float fs, float fe, float T, float dt, float phz )
/*********************************************************************
Linear - generates linear sweep
**********************************************************************
Input: 
fs	start frequency 
fe	end frequency 
T	duration in sec 
dt	sample rate in sec 
phz	initial phase

Output:
tr	array of samples
*********************************************************************
This subroutine generates a sweep with a linear frequency-time 
dependance.
*********************************************************************
References:
Any book on Vibroseis.
*********************************************************************
Author: CWP: Michel Dietrich, Tagir Galikeev		Date: 7 Oct 1994
*********************************************************************/
{
	int nt, i;
	float rate, t ;

	nt=T/dt + 1;
	rate=(fe-fs)/T;

	for (i=0; i<nt; i++) {
		t=i*dt;
		tr.data[i]=cos( TWOPI*(fs+rate/2*t)*t + phz );
	}
}

void Linear_Segment ( float *freq, float *time, int isegm, float T, 
			float dt, float phz )
/*********************************************************************
Linear_Segment - generates linear-segment sweep
**********************************************************************
Input: 
freq	array of start and end frequencies for each segment 
time	array of time lenghts for each segment
isegm	total number of segments in a sweep
T	duration in sec 
dt	sample rate in sec 
phz	initial phase

Output:
tr	array of samples
*********************************************************************
This subroutine generates a nonlinear sweep which consists of a 
number of linear segments
*********************************************************************
References:
Any book on Vibroseis.
*********************************************************************
Author: Tagir Galikeev				  Date:7 Oct 1994
*********************************************************************/
{
	int nt, m, k, i, j;
	float aa, ab, phi, phase;
	
	nt=T/dt + 1;
	m = 0;
	k = 0;
	phase = 0.;

	for (i=0; i<isegm; i++) {
		k = k + m;
		m = (int)( (float)nt / T * (float)time[i] );
		aa = TWOPI*freq[i]*dt;
		ab = PI*(freq[i+1]-freq[i])*(dt*dt) / time[i];

		for (j=0; j<m; j++) {
	  		phi = j*(ab*j+aa) + phase + phz;
			tr.data[j+k]=cos( phi );
		}

		phase = (ab*m+aa)*m + phase;

	}
}

void dB_per_octave ( float fs, float fe, float T, float dt, float swconst, 
			float phz )
/*********************************************************************
dB_per_octave - generates a decibel per octave sweep
**********************************************************************
Input: 
fs	start frequency of a sweep
fe	end frequency of a sweep
T	 duration in sec 
dt	sample rate in sec 
swconst	sweep constant (boost)
phz	initial phase

Output:
tr	array of samples
*********************************************************************
This subroutine generates a nonlinear sweep with a boost in decibels
per octave
*********************************************************************
References:
Equation is taken from PELTON manual
*********************************************************************
Author: Tagir Galikeev				  Date:7 Oct 1994
*********************************************************************/
{
	int nt, i;
	float K1, K2, t, s, swcon;

	nt=T/dt + 1;
	if( swconst == -6. ) swconst=-5.999;
	swcon=swconst/6+1;
	s=(swcon+1)/swcon;
	K1=(float)pow((double)fs,(double)swcon);
	K2=( (float)pow((double)fe,(double)swcon)-
		(float)pow((double)fs,(double)swcon) ) / T;

	for (i=0; i<nt; i++) {
		t=i*dt;
		tr.data[i]=cos(TWOPI/(s*K2)*pow((double)(K1+K2*t),
							(double)s) + phz);
	}
}

void dB_per_hertz ( float fs, float fe, float T, float dt, float swconst,
			float phz )
/*********************************************************************
dB_per_hertz - generates a decibel per hertz sweep
**********************************************************************
Input: 
fs	start frequency of a sweep
fe	end frequency of a sweep
T	 duration in sec 
dt	sample rate in sec 
swconst	sweep constant (boost)
phz	initial phase

Output:
tr	array of samples
*********************************************************************
This subroutine generates a nonlinear sweep with a boost in decibels
per hertz
*********************************************************************
References:
Equation is taken from PELTON manual
*********************************************************************
Author: Tagir Galikeev				  Date:7 Oct 1994
*********************************************************************/
{
	int nt, i;
	float K1, K2, t;

	nt=T/dt + 1;
	K1=(float)20./(swconst*log(10.));
	K2=(float)(exp(swconst*log(10.)*(fe-fs)/20.) - 1.)/T;

	for (i=0; i<nt; i++) {
		t=i*dt;
		if (swconst==0) tr.data[i]=fs;
		else tr.data[i]=cos(TWOPI*(fs*t+K1/K2*((1.+t*K2)*log(1.+t*K2)-
						(1.+t*K2))) + phz);
	}
}

void t_power ( float fs, float fe, float T, float dt, float swconst,
		float phz )
/*********************************************************************
t_power - generates a time to the power sweep
**********************************************************************
Input: 
fs	start frequency of a sweep
fe	end frequency of a sweep
T	 duration in sec 
dt	sample rate in sec 
swconst	sweep constant (boost)
phz	initial phase

Output:
tr	array of samples
*********************************************************************
This subroutine generates a nonlinear sweep with a time to
swconst power boost
*********************************************************************
References:
Equation is taken from PELTON manual
*********************************************************************
Author: Tagir Galikeev				  Date:7 Oct 1994
*********************************************************************/
{
	int nt, i;
	float t, s;

	nt=T/dt + 1;
	for (i=0; i<nt; i++) {
		t=i*dt; s=t/T;
		tr.data[i]=cos( TWOPI*t*( fs+(fe-fs)/(swconst+1.)*
				pow((double)s,(double)swconst) ) + phz );
	}
}


#define EPS     3.8090232       /* exp(-EPS*EPS) = 5e-7, "noise" level  */
				/* see sugain.c				*/
void sweep_taper ( float t1, float t2, int tap_type, float T, float dt )
/*********************************************************************
sweep taper - tapers the sweep
**********************************************************************
Input: 
t1	start taper in sec
t2	end taper in sec
tap_type  type of taper to apply: 1 linear, 2 sine, 3 cosine
T	 sweep duration in sec
dt	sample rate in sec
 
Output:
tr	array of tapered samples
*********************************************************************
This subroutine tapers a sweep mainly to reduce Gibbs phenomena.
Taper coulld be one of the specified above.
*********************************************************************
References:
Any book on Vibroseis.
*********************************************************************
Author: Tagir Galikeev				  Date:7 Oct 1994
*********************************************************************/
{
	int nt, i, nt1, nt2;
	float env=0.0, f, x;

	nt = (int)(T / dt + 1);
	nt1 = (int)(t1 / dt + 1);
	nt2 = (int)(t2 / dt + 1);
	/* apply start taper */
	if( nt1 > 1 ) {
		for (i=0; i<nt1; i++) {
	  		f = (float)i / (float)nt1;
	  		switch ((char) tap_type)	{
	  			case 1: env=f;
	  				break;
	  			case 2: env=sin(PI*f/2.);
	  				break;
	  			case 3: env=0.5*(1.0-cos(PI*f));
	  				break;
	  			case 4: x=EPS*(1-f);
	  				env=exp(-(x*x));
	  				break;
	  			case 5: x=2.0*(1-f);
	  				env=exp(-(x*x));
	  				break;
	  			default:throw cseis_geolib::csException(" taper ?!");
	  		}
	  		tr.data[i] *= env;
		}
	}
	/* apply end taper */
	if( nt2 > 1 ) {
		for (i=0; i<nt2; i++) {
	  		f = (float)i / (float)nt2;
			switch ((char) tap_type)	{
	  			case 1: env=f;
	  				break;
	  			case 2: env=sin(PI*f/2.);
	  				break;
	  			case 3: env=0.5*(1.0-cos(PI*f));
	  				break;
	  			case 4: x=EPS*(1-f);
	  				env=exp(-(x*x));
	  				break;
	  			case 5: x=2.0*(1-f);
	  				env=exp(-(x*x));
	  				break;
	  			default:throw cseis_geolib::csException(" taper ?!");
	  		}
			tr.data[nt-i]  *= env;
		}
	}
}

} // END namespace
