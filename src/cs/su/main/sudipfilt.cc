/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUDIPFILT: $Revision: 1.20 $ ; $Date: 2011/11/12 00:09:00 $		*/

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
extern "C" {
  #include <signal.h>

}

/*********************** self documentation **********************/
std::string sdoc_sudipfilt =
" 								"
" SUDIPFILT - DIP--or better--SLOPE Filter in f-k domain	"
" 								"
" sudipfilt <infile >outfile [optional parameters]		"
"								"
" Required Parameters:						"
" dt=(from header)	if not set in header then mandatory	"
" dx=(from header, d1)	if not set in header then mandatory	"
"								"
" Optional parameters:						"
" slopes=0.0		monotonically increasing slopes		"
" amps=1.0		amplitudes corresponding to slopes	"
" bias=0.0		slope made horizontal before filtering	"
"								"
" verbose=0	verbose = 1 echoes information			"
"								"
" tmpdir= 	 if non-empty, use the value as a directory path"
"		 prefix for storing temporary files; else if the"
"	         the CWP_TMPDIR environment variable is set use	"
"	         its value for the path; else use tmpfile()	"
" 								"
" Notes:							"
" d2 is an acceptable alias for dx in the getpar		"
"								"
" Slopes are defined by delta_t/delta_x, where delta		"
" means change. Units of delta_t and delta_x are the same	"
" as dt and dx. It is sometimes useful to fool the program	"
" with dx=1 dt=1, thus avoiding units and small slope values.	"
"								"
" Linear interpolation and constant extrapolation are used to	"
" determine amplitudes for slopes that are not specified.	"
" Linear moveout is removed before slope filtering to make	"
" slopes equal to bias appear horizontal.  This linear moveout	"
" is restored after filtering.  The bias parameter may be	"
" useful for spatially aliased data.  The slopes parameter is	"
" compensated for bias, so you need not change slopes when you	"
" change bias.							"
" 								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sudipfilt {


/* Credits:
 *
 *	CWP: Dave (algorithm--originally called slopef)
 *	     Jack (reformatting for SU)
 *
 * Trace header fields accessed: ns, dt, d2
 */
/**************** end self doc ***********************************/


#define PFA_MAX	720720	/* Largest allowed nfft	          */

/* prototypes */
void slopefilter (int nslopes, float slopes[], float amps[], float bias,
	int nt, float dt, int ntr, float dx, FILE *datafp);
static void closefiles(void);

/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/

segy tr;

void* main_sudipfilt( void* args )
{
	int nt;		/* number of time samples			*/
	float dt;	/* time sampling interval			*/
	int ntr;	/* number of traces				*/
	float dx;	/* trace spacing (spatial sampling interval)	*/
	int nslopes;	/* number of slopes specified			*/
	float *slopes;	/* slopes at which amplitudes are specified	*/
	int namps;	/* number of amplitudes specified		*/
	float *amps;	/* amplitudes corresponding to slopes		*/
	float bias;	/* slope bias					*/
	int verbose;	/* flag for echoing info			*/
	char *tmpdir;	/* directory path for tmp files			*/
	cwp_Bool istmpdir=cwp_false;/* true for user-given path		*/


	/* Hook up getpar to handle the parameters */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sudipfilt );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get info from first trace */ 
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;


	/* Get parameters */
	if (!parObj.getparint("verbose", &verbose))	verbose = 0;

	/* Look for user-supplied tmpdir */
	if (!parObj.getparstring("tmpdir",&tmpdir) &&
	    !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
	if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
		throw cseis_geolib::csException("you can't write in %s (or it doesn't exist)", tmpdir);

	if (!parObj.getparfloat("dt", &dt))	dt = ((double) tr.dt)/1000000.0;
	if (!dt) throw cseis_geolib::csException("dt field is zero and not parObj.getparred");
	if (!parObj.getparfloat("dx", &dx) && !getparfloat("d2", &dx))	dx = tr.d2;
	if (!dx) throw cseis_geolib::csException("d2 field is zero and dx not parObj.getparred");

	nslopes = parObj.countparval("slopes");
	if (nslopes) {
		slopes = alloc1float(nslopes);
		parObj.getparfloat("slopes", slopes);
	} else {
		nslopes = 1;
		slopes = alloc1float(nslopes);
		slopes[0] = 0.0;
	}

	namps = parObj.countparval("amps");
	if (namps) {
		amps = alloc1float(namps);
		parObj.getparfloat("amps", amps);
	} else {
		namps = 1;
		amps = alloc1float(namps);
		amps[0] = 1.0;
		warn("no amps given--doing no-op");
	}

	if (!parObj.getparfloat("bias", &bias)) bias = 0.0;


	/* Check parameters */
	if (nslopes != namps)
		throw cseis_geolib::csException("number of slopes (%d) must equal number of amps(%d)",
			nslopes, namps);
	{ register int i;
	  for (i=1; i<nslopes; ++i)
		if (slopes[i] <= slopes[i-1])
			throw cseis_geolib::csException("slopes must be monotonically increasing");
	}


	/* Store traces and headers in tmpfile while getting a count */
	if (STREQ(tmpdir,"")) {
		tracefp = etmpfile();
		headerfp = etmpfile();
		if (verbose) warn("using tmpfile() call");
	} else { /* user-supplied tmpdir */
		char directory[BUFSIZ];
		strcpy(directory, tmpdir);
		strcpy(tracefile, temporary_filename(directory));
		strcpy(headerfile, temporary_filename(directory));
		/* Trap signals so can remove temp files */
		signal(SIGINT,  (void (*) (int)) closefiles);
		signal(SIGQUIT, (void (*) (int)) closefiles);
		signal(SIGHUP,  (void (*) (int)) closefiles);
		signal(SIGTERM, (void (*) (int)) closefiles);
		tracefp = efopen(tracefile, "w+");
		headerfp = efopen(headerfile, "w+");
      		istmpdir=cwp_true;		
		if (verbose) warn("putting temporary files in %s", directory);
	}

	ntr = 0;
	do { 
		++ntr;
		efwrite(&tr, 1, HDRBYTES, headerfp);
		efwrite(tr.data, FSIZE, nt, tracefp);
	} while (cs2su->getTrace(&tr));
	

	/* Apply slope filter */
	slopefilter(nslopes,slopes,amps,bias,nt,dt,ntr,dx,tracefp);


	/* Output filtered traces */
	erewind(headerfp);
	erewind(tracefp);
	{ register int itr;
	  for (itr = 0; itr < ntr; ++itr) {
		efread(&tr, 1, HDRBYTES, headerfp);
		efread(tr.data, FSIZE, nt, tracefp);
		su2cs->putTrace(&tr);
	  }
	}


	/* Clean up */
	efclose(headerfp);
	if (istmpdir) eremove(headerfile);
	efclose(tracefp);
	if (istmpdir) eremove(tracefile);
	free1float(slopes);
	free1float(amps);

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



void slopefilter (int nslopes, float slopes[], float amps[], float bias,
	int nt, float dt, int nx, float dx, FILE *tracefp)
/******************************************************************************
apply slope filter in frequency-wavenumber domain
*******************************************************************************
Input:
nslopes		number of slopes (and amplitudes) specified
slopes		slopes at which amplitudes are specified (see notes below)
amps		amplitudes corresponding to slopes (see notes below)
bias		linear moveout slope before and after filtering
nt		number of time samples
dt		time sampling interval
nx		number of traces
dx		trace space (spatial sampling interval)
tracefp		file pointer to data to be filtered

Output:
tracefp		file pointer to filtered data
*******************************************************************************
Notes:
Linear interpolation and constant extrapolation are used to
determine amplitudes for slopes that are not specified.
******************************************************************************/
{
	int ntfft;		/* nt after padding for FFT */
	int nxfft;		/* nx after padding for FFT */
	float sfft;		/* scale factor for FFT */
	int nw;			/* number of frequencies */
	float dw;		/* frequency sampling interval */
	float fw;		/* first frequency */
	int nk;			/* number of wavenumbers */
	float dk;		/* wavenumber sampling interval */
	float w,k;		/* frequency and wavenumber */
	int it,ix,iw,ik;	/* sample indices */
	float slope,amp;	/* slope and amplitude for particular w,k */
	complex **cpfft;	/* complex FFT workspace */
	float **pfft;		/* float FFT workspace */
	float phase;		/* phase shift for bias */
	complex cshift;		/* complex phase shifter for bias */

	/* determine lengths and scale factors for prime-factor FFTs */
	ntfft = npfar(nt);
	nxfft = npfa(nx);
	sfft = 1.0/(ntfft*nxfft);
	
	/* determine frequency and wavenumber sampling */
	nw = ntfft/2+1;
	dw = 2.0*PI/(ntfft*dt);
	fw = 0.000001*dw; /* non-zero to avoid divide by zero w */
	nk = nxfft;
	dk = 2.0*PI/(nxfft*dx);

	/* allocate real and complex workspace for FFTs */
	cpfft = alloc2complex(nw,nk);
	pfft = alloc2float(ntfft,nxfft);

	/* copy data from input to FFT array and pad with zeros */
	rewind(tracefp);
	for (ix=0; ix<nx; ix++) {
		efread(pfft[ix], FSIZE, nt, tracefp);
		for (it=nt; it<ntfft; it++)
			pfft[ix][it] = 0.0;
	}
	for (ix=nx; ix<nxfft; ix++)
		for (it=0; it<ntfft; it++)
			pfft[ix][it] = 0.0;
	
	/* Fourier transform t to w */
	pfa2rc(1,1,ntfft,nx,pfft[0],cpfft[0]);
	
	/* do linear moveout bias via phase shift */
	for (ix=0; ix<nx; ix++) {
		for (iw=0,w=0.0; iw<nw; iw++,w+=dw) {
			phase = -ix*dx*w*bias;
			cshift = cmplx(cos(phase),sin(phase));
			cpfft[ix][iw] = cmul(cpfft[ix][iw],cshift);
		}
	}
	
	/* Fourier transform x to k */
	pfa2cc(-1,2,nw,nxfft,cpfft[0]);
	
	/* loop over wavenumbers */
	for (ik=0; ik<nk; ik++) {
	
		/* determine wavenumber */
		k = (ik<=nk/2) ? ik*dk : (ik-nk)*dk;
		
		/* loop over frequencies */
		for (iw=0,w=fw; iw<nw; iw++,w+=dw) {
		
			/* determine biased slope */
			slope = k/w+bias;
			
			/* linearly interpolate to find amplitude */
			intlin(nslopes,slopes,amps,amps[0],amps[nslopes-1],
				1,&slope,&amp);
			
			/* include fft scaling */
			amp *= sfft;
			
			/* filter real and imaginary parts */
			cpfft[ik][iw].r *= amp;
			cpfft[ik][iw].i *= amp;
		}
	}

	/* Fourier transform k to x */
	pfa2cc(1,2,nw,nxfft,cpfft[0]);

	/* undo linear moveout bias via phase shift */
	for (ix=0; ix<nx; ix++) {
		for (iw=0,w=0.0; iw<nw; iw++,w+=dw) {
			phase = ix*dx*w*bias;
			cshift = cmplx(cos(phase),sin(phase));
			cpfft[ix][iw] = cmul(cpfft[ix][iw],cshift);
		}
	}

	/* Fourier transform w to t */
	pfa2cr(-1,1,ntfft,nx,cpfft[0],pfft[0]);
	
	/* copy filtered data from FFT array to output */
	rewind(tracefp);
	for (ix=0; ix<nx; ix++)
		efwrite(pfft[ix], FSIZE, nt, tracefp);

	/* free workspace */
	free2complex(cpfft);
	free2float(pfft);
}

/* for graceful interrupt termination */
static void closefiles(void)
{
	efclose(headerfp);
	efclose(tracefp);
	eremove(headerfile);
	eremove(tracefile);
	exit(EXIT_FAILURE);
}

} // END namespace
