/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUKFILTER: $Revision: 1.7 $ ; $Date: 2011/11/12 00:09:00 $	*/

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

/*********************** self documentation **********************/
std::string sdoc_sukfilter =
" 									"
" SUKFILTER - radially symmetric K-domain, sin^2-tapered, polygonal	"
"		  filter						"
" 									"
"     sukfilter <infile >outfile [optional parameters]			"
"									"
" Optional parameters:							"
" k=val1,val2,...	array of K filter wavenumbers			"
" amps=a1,a2,...	array of K filter amplitudes			"
" d1=tr.d1 or 1.0	sampling interval in first (fast) dimension	"
" d2=tr.d1 or 1.0	sampling interval in second (slow) dimension	"
"									"
" Defaults:								"
" k=.10*(nyq),.15*(nyq),.45*(nyq),.50*(nyq)				"
" amps=0.,1.,...,1.,0.  trapezoid-like bandpass filter			"
"									"
" The nyquist wavenumbers, nyq=sqrt(nyq1^2 + nyq2^2) is  computed	"
" internally.								"
"									"
" Notes:								"
" The filter is assumed to be symmetric, to yield real output.		"
"									"
" Because the data are assumed to be purely spatial (i.e. non-seismic), "
" the data are assumed to have trace id (30), corresponding to (z,x) data"
"									"
" The relation: w = 2 pi F is well known for frequency, but there	"
" doesn't seem to be a commonly used letter corresponding to F for the	"
" spatial conjugate transform variables.  We use K1 and K2 for this.	"
" More specifically we assume a phase:					"
"		-i(k1 x1 + k2 x2) = -2 pi i(K1 x1 + K2 x2).		"
" and K1, K2 define our respective wavenumbers.				"
" 									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sukfilter {


/* Credits:
 *     CWP: John Stockwell, June 1997.
 *
 * Trace header fields accessed: ns, d1, d2
 */
/**************** end self doc ***********************************/


/* definitions */
#define PFA_MAX	720720	/* Largest allowed nfft		  */
#define LOOKFAC 2	/* look factor			  */
#define FRAC0   0.10    /* Ratio of default k1 to Nyquist */
#define FRAC1   0.15    /* Ratio of default k2 to Nyquist */
#define FRAC2   0.45    /* Ratio of default k3 to Nyquist */
#define FRAC3   0.50    /* Ratio of default k4 to Nyquist */

#define N_PHI	3	/* arbitrary number of phi values */

/* prototype of function used internally */
void polygonalFilter(float *f, float *amps, int npoly,
				int nfft, float dt, float *filter);
segy tr;

void* main_sukfilter( void* args )
{
	int nx1,nx2;		/* numbers of samples			*/
	float dx1,dx2,dx;	/* sampling intervals			*/
	int verbose;		/* verbose flag				*/
       
        float nyq;		/* K Nyquist wavenumber			*/

	int npoly;		/* number of points defining k filter	*/
	int namps;		/* number of amps defining k filter	*/

        float *k;		/* wavenumber values defining k filter	*/
	float **kphi;		/* k,phi array				*/
	int nphi=1;		/* number of phi values for polartorect */
	int iphi;		/* phi values				*/
        float *amps;		/* amplitude values defining k filter	*/

	float *kfilt;		/* k wavenumber filter			*/
	float **karray;		/* array of filter values 		*/
	float **kfilter;	/* wavenumber filter 			*/

	int ix1,ix2;		/* sample indices			*/
	int nx1fft,nx2fft,nxfft;/* dimensions after padding for FFT	*/
	int nK1,nK2,nK;		/* transform (output) dimensions	*/

        int ik;			/* k counter				*/
        int ik1;		/* k1 counter				*/
        int ik2;		/* k2 counter				*/
        int iamps;		/* amplitude counter			*/
        int icount;		/* zero counter				*/

	register complex **ct;	/* complex FFT workspace		*/
	register float **rt;	/* float FFT workspace			*/
	FILE *tracefp;		/* temp file to hold traces		*/
	FILE *hfp;		/* temp file to hold trace headers	*/

	/* Hook up getpar to handle the parameters */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sukfilter );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get info from first trace */ 
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	if (tr.trid != TRID_DEPTH)
				warn("tr.trid = %d",tr.trid);
	nx1 = tr.ns;

	/* get sampling intervals */
	if (!parObj.getparfloat("d1", &dx1)) {
		if (tr.d1) { /* is dt field set? */
			dx1 = (float) tr.d1;
		} else { /* d1 not set, assume 1.0 */
			dx1 = 1.0;
			warn("tr.d1 not set, assuming d1=1.0");
		}
	}
	if (!parObj.getparfloat("d2",&dx2)) {
		if (tr.d2) { /* is d2 field set? */
			dx2 = tr.d2;
		} else {
			dx2 = 1.0;
			warn("tr.d2 not set, assuming d2=1.0");
		}
	}
	dx = NINT(sqrt(dx1*dx1 + dx2*dx2));

	if (!parObj.getparint("verbose",&verbose))	verbose = 1;


	/* Compute Nyquist wavenumber */
	nyq = 0.5/dx;

	/* Store traces in tmpfile while getting a count */
	tracefp = etmpfile();
	hfp = etmpfile();
	nx2 = 0;
	do { 
		++nx2;
		efwrite(&tr, HDRBYTES, 1, hfp);
		efwrite(tr.data, FSIZE, nx1, tracefp);
	} while (cs2su->getTrace(&tr));


	/* Determine lengths for prime-factor FFTs */
	nx1fft = npfaro(nx1, LOOKFAC*nx1);
	nx2fft = npfa(nx2);
	if (nx1fft >= SU_NFLTS || nx1fft >= PFA_MAX)
		throw cseis_geolib::csException("Padded nx1=%d--too big",nx1fft);
	if (nx2fft >= PFA_MAX)
		throw cseis_geolib::csException("Padded nx2=%d--too big",nx2fft);

	/* Determine number of wavenumbers in K1 and K2 */
	nK1 = nx1fft/2 + 1;
	nK2 = nx2fft/2 + 1;
	nK = NINT(sqrt(nK1*nK1 + nK2*nK2));
	nxfft = 2*( nK - 1); /* faked for polygonal filter in k */

        /* Get wavenumbers that define the k filter */
        if ((npoly = parObj.countparval("k"))!=0) {
                k = ealloc1float(npoly);
                parObj.getparfloat("k",k);
        } else {
                npoly = 4;
                k = ealloc1float(npoly);

                k[0] = FRAC0 * nyq;
                k[1] = FRAC1 * nyq;
                k[2] = FRAC2 * nyq;
                k[3] = FRAC3 * nyq;
        }

	/* Check k values */
	if(npoly < 2) warn("Only %d value defining filter",npoly);
        for(ik=0; ik < npoly-1; ik++)
		if(k[ik] < 0.0 || k[ik] > k[ik+1])
                                throw cseis_geolib::csException("Bad filter parameters");

	/* Get k filter amplitude values */
        if ((namps = parObj.countparval("amps"))!=0) {
                amps = ealloc1float(namps);
                parObj.getparfloat("amps",amps);
        } else {
                namps = npoly;
                amps = ealloc1float(namps);

		/* default is a trapezoidal bandpass filter */
		for(iamps=1; iamps<namps-1; iamps++) amps[iamps]=1.;

		amps[0]=0.;
		amps[namps-1]=0.;
        }
	if (!(namps==npoly)) 
		throw cseis_geolib::csException("number of k values must = number of amps values");

        /* Check amps values */
        for(iamps = 0, icount=0; iamps < namps ; iamps++) {
		icount+=amps[iamps];
                if( amps[iamps] < 0.) throw cseis_geolib::csException("amps values must be positive");
        }
        if (icount==0) throw cseis_geolib::csException("All amps values are zero");
        for(iamps = 0, icount=0; iamps < namps-1 ; ++iamps) {
			if(!(amps[iamps]==amps[iamps+1])) ++icount;
	}
        if (icount==0) warn("All amps values are the same");

	/* Allocate space */
	rt = ealloc2float(nx1fft, nx2fft);
	ct = ealloc2complex(nK1,nx2fft);
	karray = ealloc2float(nx1fft,nx2fft);
	kfilter = ealloc2float(nx1fft,nx2fft);
	kfilt = ealloc1float(nK);
	kphi = ealloc2float(1,nK);

	/* Zero all arrays */
	memset((void *) rt[0], 0, nx1fft*nx2fft*FSIZE);
	memset((void *) kfilter[0], 0, nx1fft*nx2fft*FSIZE);
	memset((void *) ct[0], 0, nK1*nx2fft*sizeof(complex));
	memset((void *) kfilt, 0, nK*FSIZE);
	memset((void *) kphi[0], 0, nK*FSIZE);

	/* Build Filter */
	polygonalFilter(k,amps,npoly,nxfft,dx,kfilt);

	/* Build k,phi filter */
	for (iphi=0; iphi<1 ; ++iphi)
		for (ik=0; ik<nK; ++ik)
			kphi[ik][iphi] = kfilt[ik]*kfilt[ik];

	/* convert polar to rectangular coordinates */
	polartorect(nphi,0,0,nK,1,0,
			kphi,nx1fft,1,-nK1,nx2fft,
			1,-nK2,karray);

	/* positive k1, positive k2 */
	for (ik2=0; ik2<nK2; ++ik2) 
		for (ik1=0; ik1<nK1; ++ik1) 
			kfilter[ik2][ik1]=karray[nK2-1-ik2][nK1-1-ik1];

	/* positive k1, negative k2 */
	for (ik2=nK2; ik2<nx2fft; ++ik2) 
		for (ik1=0; ik1<nK1; ++ik1) 
			kfilter[ik2][ik1]=karray[ik2-nK2][nK1-1-ik1];

	/* Load traces into fft arrays and close tmpfile */
	rewind(tracefp);
	for (ix2=0; ix2<nx2; ++ix2)
		efread(rt[ix2], FSIZE, nx1, tracefp);

	efclose(tracefp);
	
	/* Fourier transform dimension 1 */
	pfa2rc(-1,1,nx1fft,nx2,rt[0],ct[0]);

	/* Fourier transform dimension 2 */
	pfa2cc(-1,2,nK1,nx2fft,ct[0]);

	/* Apply filter */
	for (ik2=0; ik2 < nx2fft; ++ik2)
		for (ik1=0; ik1 < nK1; ++ik1)
			ct[ik2][ik1] = crmul(ct[ik2][ik1], kfilter[ik2][ik1]);


	/* Inverse Fourier transform dimension 2 */
	pfa2cc(1,2,nK1,nx2fft,ct[0]);
	
	/* Inverse Fourier transform dimension 1 */
	pfa2cr(1,1,nx1fft,nx2,ct[0],rt[0]);

	erewind(hfp);
	/* Output filtered traces */
	for (ix2=0; ix2 < nx2; ++ix2) { 

		efread(&tr, HDRBYTES, 1, hfp);
		for (ix1=0; ix1 <nx1 ; ++ix1) 
			tr.data[ix1] = rt[ix2][ix1];

		su2cs->putTrace(&tr);
	}

	efclose(hfp);

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

void polygonalFilter(float *f, float *amps, int npoly, 
				int nfft, float dt, float *filter)
/*************************************************************************
polygonalFilter -- polygonal filter with sin tapering
**************************************************************************
Input:
f		array[npoly] of frequencies defining the filter
amps		array[npoly] of amplitude values
npoly		size of input f and amps arrays
dt		time sampling interval
nfft		number of points in the fft
pow		power flag, use pow=1 for 1D usage, =2 for 2D

Output:
filter		array[nfft] filter values
**************************************************************************
Notes: Filter is to be applied in the frequency domain
Also, this is not exactly the same filter that is used in "sufilter".
That filter is sin^2. This done because the filter is squared in
its application to make the amplitudes come out right.
**************************************************************************
Author:  CWP: John Stockwell   1992
*************************************************************************/
#define PIBY2   1.57079632679490
{
        int *intfr;             /* .... integerizations of f		*/
        int icount,ifs;		/* loop counting variables              */
	int taper=0;		/* flag counter				*/
        int nf;                 /* number of frequencies (incl Nyq)     */
        int nfm1;               /* nf-1                                 */
        float onfft;            /* reciprocal of nfft                   */
        float df;               /* frequency spacing (from dt)          */

	intfr=alloc1int(npoly);

        nf = nfft/2 + 1;
        nfm1 = nf - 1;
        onfft = 1.0 / nfft;

        /* Compute array of integerized frequencies that define the filter*/
        df = onfft / dt;
        for(ifs=0; ifs < npoly ; ++ifs) {
                intfr[ifs] = NINT(f[ifs]/df);
                if (intfr[ifs] > nfm1) intfr[ifs] = nfm1;
        }

	/* Build filter, with scale, and taper specified by amps[] values*/
	/* Do low frequency end first*/
	for(icount=0; icount < intfr[0] ; ++icount) 
		filter[icount] = amps[0] * onfft;

	/* now do the middle frequencies */
	for(ifs=0 ; ifs<npoly-1 ; ++ifs){
	   if(amps[ifs] < amps[ifs+1]) {	
		++taper;
		for(icount=intfr[ifs]; icount<=intfr[ifs+1]; ++icount) {
		    float c = PIBY2 / (intfr[ifs+1] - intfr[ifs] + 2);
		    float s = sin(c*(icount - intfr[ifs] + 1));
		    float adiff = amps[ifs+1] - amps[ifs];
		    filter[icount] = (amps[ifs] + adiff*s) * onfft;
		}
	   } else if (amps[ifs] > amps[ifs+1]) {	
		++taper;
		for(icount=intfr[ifs]; icount<=intfr[ifs+1]; ++icount) {
			   float c = PIBY2 / (intfr[ifs+1] - intfr[ifs] + 2);
                	   float s = sin(c*(intfr[ifs+1] - icount + 1));
			   float adiff = amps[ifs] - amps[ifs+1];
                	   filter[icount] = (amps[ifs+1] + adiff*s) * onfft;
		  }
	   } else 
		if(!(taper)){
		for(icount=intfr[ifs]; icount <= intfr[ifs+1]; ++icount)
		   	   filter[icount] = amps[ifs] * onfft;
		} else {
		for(icount=intfr[ifs]+1; icount <= intfr[ifs+1]; ++icount)
		   	   filter[icount] = amps[ifs] * onfft;
		}
	}

	/* finally do the high frequency end */
	for(icount=intfr[npoly-1]+1; icount<nf; ++icount){
		filter[icount] = amps[npoly-1] * onfft;
	}

}

} // END namespace
