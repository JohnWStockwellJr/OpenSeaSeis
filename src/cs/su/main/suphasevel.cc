/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUPHASEVEL: $Revision: 1.3 $ ; $Date: 2011/11/16 23:35:04 $		*/

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
extern "C" {
  #include <signal.h>

}

/*********************** self documentation **********************/
std::string sdoc_suphasevel =
" 								"
" SUPHASEVEL - Multi-mode PHASE VELocity dispersion map computed"
"              from shot record(s)				"
" 								"
" suphasevel <infile >outfile [optional parameters]		"
"								"
" Optional parameters:						"
" fv=330	minimum phase velocity (m/s)			"
" nv=100	number of phase velocities			"
" dv=25		phase velocity step (m/s)			"
" fmax=50	maximum frequency to process (Hz)		"
"		=0 process to nyquist				"
" norm=0	do not normalize by amplitude spectrum		"
"		=1 normalize by amplitude spectrum		"
" verbose=0	verbose = 1 echoes information			"
" 								"
" Notes:  Offsets read from headers.			 	"
"  1. output is complex frequency data				"
"  2. offset header word must be set (signed offset is ok)	"
"  3. norm=1 tends to blow up aliasing and other artifacts	"
"  4. For correct suifft later, use fmax=0			"
"  5. For later processing outtrace.dt=domega			"
"  6. works for 2D or 3D shots in any offset order		"
"								"
" Example:							"
"   suphasevel < shotrecord.su | suamp | suximage 	  	"
"								"
" Reference: Park, Miller, and Xia (1998, SEG Abstracts)	"
"								"
" Trace header fields accessed: dt, offset, ns			"
" Trace header fields modified: nx,dt,trid,d1,f1,d2,f2,tracl	"
" 								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suphasevel {


/* Credits:
 *
 *	UHouston: Chris Liner June2008 (cloned from suspecfk)
 *
 *  This code implements the following integral transform
 *             _
 *            /
 *  u(w,v) = / k(w,x,v) u(w,x) dx
 *         _/
 *  where
 *	u(w,v) is the phase velocity dispersion image
 *	k(w,x,v) is the transform kernel.... exp(-i w x / v)
 *	u(w,x) = FT[u(t,x)] is the input shot record(s) 
 *
 * Trace header fields accessed: dt, offset, ns
 * Trace header fields modified: nx,dt,trid,d1,f1,d2,f2,tracl
 */
/**************** end self doc ***********************************/

#define	I		cmplx(0.0, 1.0)
#define PFA_MAX	720720	/* Largest allowed nfft	          	*/
#define MAX_OFFS 10000 	/* Largest allowable offset		*/

/* Prototype */
static void closefiles(void);

/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
FILE *tracefp;		/* fp for trace storage file		*/

segy intrace, outtrace;

void* main_suphasevel( void* args )
{
	int nt,nx;		/* numbers of samples			*/
	float dt;		/* sampling intervals			*/
	int it,ix;		/* sample indices			*/
	int ntfft;		/* dimensions after padding for FFT	*/
	int nF;			/* transform (output) dimensions	*/
	int iF;			/* transform sample indices		*/

	register complex **ct=NULL;	/* complex FFT workspace	*/
	register float **rt=NULL;	/* float FFT workspace		*/

	int verbose;		/* flag for echoing information		*/

	char *tmpdir=NULL;	/* directory path for tmp files		*/
	cwp_Bool istmpdir=cwp_false;/* true for user-given path		*/

	float v,fv,dv;		/* phase velocity, first, step	 	*/
	float amp,oamp;		/* temp vars for amplitude spectrum	*/
	int nv,iv;		/* number of phase vels, counter	*/
	float x;		/* offset  				*/
	float omega;		/* circular frequency			*/
	float domega;		/* circular frequency spacing (from dt)	*/
	float onfft;		/* 1 / nfft				*/
	float phi;		/* omega/phase_velocity			*/
	complex *cDisp=NULL;	/* temp array for complex dispersion	*/
	float arg;		/* temp var for phase calculation	*/
	complex cExp;		/* temp vars for phase calculation	*/
	float *offs=NULL;	/* input data offsets			*/
	float fmax;		/* max freq to proc (Hz)    		*/

	int out;		/* output real or abs v(f) spectrum	*/
	int norm;		/* normalization flag			*/

	float xmax;		/* maximum abs(offset) of input		*/
	float twopi, f;		/* constant and frequency (Hz)		*/



	/* Hook up getpar to handle the parameters */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suphasevel );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get info from first trace */ 
	if (!cs2su->getTrace(&intrace))  throw cseis_geolib::csException("can't get first trace");
	nt = intrace.ns;

	/* dt is used only to set output header value d1 */
	if (!parObj.getparfloat("dt", &dt)) {
		if (intrace.dt) { /* is dt field set? */
			dt = ((double) intrace.dt)/ 1000000.0;
		} else { /* dt not set, exit */
			throw cseis_geolib::csException("tr.dt not set, stop.");
		}
	}
	warn("dt=%f",dt);
	if (!parObj.getparfloat("fv",&fv))	fv   = 330;
	if (!parObj.getparfloat("dv",&dv))     dv   = 25;
	if (!parObj.getparint("nv",&nv))       nv   = 100;
	if (!parObj.getparint("out",&out))     out  = 0;
	if (!parObj.getparint("norm",&norm))   norm = 0;
	if (!parObj.getparfloat("fmax",&fmax)) fmax = 50;

	if (!parObj.getparint("verbose", &verbose))	verbose = 0;

	/* Look for user-supplied tmpdir */
	if (!parObj.getparstring("tmpdir",&tmpdir) &&
	    !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
	if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
		throw cseis_geolib::csException("you can't write in %s (or it doesn't exist)", tmpdir);


        parObj.checkpars();

	/* Set up tmpfile */
	if (STREQ(tmpdir,"")) {
		tracefp = etmpfile();
		if (verbose) warn("using tmpfile() call");
	} else { /* user-supplied tmpdir */
		char directory[BUFSIZ];
		strcpy(directory, tmpdir);
		strcpy(tracefile, temporary_filename(directory));
		/* Trap signals so can remove temp files */
		signal(SIGINT,  (void (*) (int)) closefiles);
		signal(SIGQUIT, (void (*) (int)) closefiles);
		signal(SIGHUP,  (void (*) (int)) closefiles);
		signal(SIGTERM, (void (*) (int)) closefiles);
		tracefp = efopen(tracefile, "w+");
      		istmpdir=cwp_true;		
		if (verbose) warn("putting temporary files in %s", directory);
	}
	
	/* we have to allocate offs(nx) before we know nx */
	offs = alloc1float(MAX_OFFS);	

	ix = 0;
	nx = 0;
	xmax = 0.0;
	
	/* get nx and max abs(offset) */
	do { 
		++nx;
		efwrite(intrace.data, FSIZE, nt, tracefp);
		offs[ix] = intrace.offset;
		if ( abs(intrace.offset) > xmax ) xmax = abs(intrace.offset);
		++ix;
	} while (cs2su->getTrace(&intrace));
	
	/* confirm that offsets are set */
	if ( xmax == 0.0 ) throw cseis_geolib::csException("tr.offset not set, stop.");


	/* Determine lengths for prime-factor FFTs */
	ntfft = npfar(nt);
	if (ntfft >= SU_NFLTS || ntfft >= PFA_MAX)
			throw cseis_geolib::csException("Padded nt=%d--too big",ntfft);

	/* Determine complex transform sizes */
	nF = ntfft/2+1;  /* must be this nF for fft */
        onfft = 1.0 / ntfft;
	twopi = 2.0 * PI;
	domega = twopi * onfft / dt;

	/* Allocate space */
	ct = alloc2complex(nF,nx);
	rt = alloc2float(ntfft,nx);

	/* Load traces into fft arrays and close tmpfile */
	erewind(tracefp);
	for (ix=0; ix<nx; ++ix) {

		efread(rt[ix], FSIZE, nt, tracefp);

		/* pad dimension 1 with zeros */
		for (it=nt; it<ntfft; ++it)  rt[ix][it] = 0.0;
	}
	efclose(tracefp);
	
	/* Fourier transform dimension 1 */
	pfa2rc(1,1,ntfft,nx,rt[0],ct[0]);

	/* set nF for processing */
	if (fmax == 0) { 	
		/* process to nyquist */
		nF = ntfft/2+1;
	} else {
		/* process to given fmax */
		nF = (int) (0.5 * twopi * fmax / domega);
	}
	
	/* data now in (w,x) domain 
	   allocate arrays  */
	cDisp = alloc1complex(nF);	
	
	/* if requested, normalize by amplitude spectrum 
	    (normalizing by amplitude blows up aliasing and other artifacts) */			
	if (norm == 1) {
		for (iF=0; iF<nF; ++iF)  {
			/* calc this frequency */
			omega = iF * domega;
			f = omega / twopi;
			/* loop over traces */
			for (ix=0; ix<nx; ++ix) {
				/* calc amplitude at this (f,x) location */
				amp = rcabs(ct[ix][iF]);
				oamp = 1.0/amp;
				/* scale field by amp spectrum */
				ct[ix][iF] = crmul(ct[ix][iF],oamp);
			}
		}
	}
	
	/* set global output trace headers */
	outtrace.ns = 2 * nF;
	outtrace.dt = dt*1000000.;  
	outtrace.trid = FUNPACKNYQ;
	outtrace.d1 = 2.0 / (ntfft * dt); /* Hz */
	outtrace.f1 = 0;
	outtrace.d2 = dv;
	outtrace.f2 = fv;

	/* loop over phase velocities */
	for (iv=0; iv<nv; ++iv) {

		/* this velocity */
		v = fv + iv*dv;
	
		/* loop over frequencies */
		for (iF=0; iF<nF; ++iF)  {

			/* this frequency and phase */
			omega = iF * domega;
			f = omega / twopi;
			phi = omega / v;

			/* initialize */
			cDisp[iF] = cmplx(0.0,0.0);		

			/* sum over abs offset (this is ok for 3D, too) */
			for (ix=0; ix<nx; ++ix) {

				/* get this x */
				x = abs(offs[ix]);

				/* target phase */
				arg = - phi * x;
				cExp = cwp_cexp(crmul(cmplx(0.0,1.0), arg));
				
				/* phase vel profile for this frequency */				 
				cDisp[iF] = cadd(cDisp[iF],cmul(ct[ix][iF],cExp));
			}
			
		}
		
		/* set trace counter */
		outtrace.tracl = iv + 1;
			
		/* copy results to output trace 
		   interleaved format like sufft.c */
		for (iF = 0; iF < nF; ++iF) {
			outtrace.data[2*iF]   = cDisp[iF].r;
			outtrace.data[2*iF+1] = cDisp[iF].i;
		}
		
		/* output freqs at this vel */
		su2cs->putTrace(&outtrace);

	}  /* next frequency */
	
	
	/* Clean up */
	if (istmpdir) eremove(tracefile);
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

/* for graceful interrupt termination */
static void closefiles(void)
{
	efclose(tracefp);
	eremove(tracefile);
	exit(EXIT_FAILURE);
}

} // END namespace
