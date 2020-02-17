/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUKFRAC: $Revision: 1.15 $ ; $Date: 2011/12/21 23:25:42 $	*/

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
std::string sdoc_sukfrac =
" 									"
" SUKFRAC - apply FRACtional powers of i|k| to data, with phase shift "
" 									"
"     sukfrac <infile >outfile [optional parameters]			"
"									"
" Optional parameters:							"
"  power=0		exponent of (i*sqrt(k1^2 + k2^2))^power		"
"			=0 ===> phase shift only			"
"			>0 ===> differentiation				"
"			<0 ===> integration				"
"  sign=1		sign on transform exponent       		"
"  d1=1.0		x1 sampling interval				"
"  d2=1.0		x2 sampling interval				"
"  phasefac=0		phase shift by phase=phasefac*PI		"
"  ...directional derivative, active only if angle= is set ....		"
"  angle=		if set applies operator directionally in the	"
"			direction specified by the value of angle,	"
"			-360.0 <= angle < 359.99999 degress		"
"									"
" Notes:								"
" The relation: w = 2 pi F is well known for frequency, but there	"
" doesn't seem to be a commonly used letter corresponding to F for the	"
" spatial conjugate transform variables.  We use K1 and K2 for this.	"
" More specifically we assume a phase:					"
"		-i(k1 x1 + k2 x2) = -2 pi i(K1 x1 + K2 x2).		"
" and K1, K2 define our respective wavenumbers.				"
" 									"
" Algorithms 								"
" 	g(x1,x2)=Re[2DINVFFT{ ( (sign) i |k|)^power 2DFFT(f)}e^i(phase)]"
" where: 								"
"       |k| = sqrt[ (k1)^2 + (k2)^2 ] 					"
" or when angle= is set 						"
"       |k| = sqrt[ (k1  cos(angle) )^2 + (k2 sin(angle) )^2 ] 		"
" 									"
" In the default mode a factor of (i|k|)^(power) is applied in the 	"
" transform domain. For time data the time axis direction is taken to 	"
" be the k1-direction.  The effect of this filter is to differentiate   "
" the input in the normal direction to any curvilinear features in	"
" in the data, and thus be a non-directional-specific edge enhancer.	"
" 									"
" If angle= is set, then the intended effect is a derivative in the 	"
" direction specified by the angle, with the k1-direction being angle=0,"
" corresponding to curves whose normals lie in the x1-direction. 	"
" 									"
" Caveat:								"
" Large amplitude errors will result of the data set has too few points."
" 									"
" Examples:								"
" Edge sharpening: 							"
" Laplacean : 								"
"    sukfrac < image_data  power=2 phasefac=-1 | ... 			"
" 									"
" Image enhancement:							"
"   Derivative filter:							"
"    sukfrac < image_data  power=1 phasefac=-.5 | ... 			"
" 									"
" Image enhancement:							"
"   Half derivative (this one is the best for photographs): 		"
"    sukfrac < image_data  power=.5 phasefac=-.25 | ... 		"
" 									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sukfrac {


/* Credits:
 *     CWP: John Stockwell, June 1997, based on sufrac.
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

#define I	cmplx(0.0,1.0)
#define PIBY2	(0.5 * PI)
#define TWOPI	(2.0 * PI)
#define SQRT2	(sqrt(2.0))

segy tr;

void* main_sukfrac( void* args )
{
	int nx1,nx2;		/* numbers of samples			*/
	float dx1,dx2;		/* sampling intervals			*/
	int verbose;		/* verbose flag				*/

	float k=0.0;		/* magnitude of wavenumber		*/
	float amp;		/* amplitude associated with power	*/
	complex **kfilter;	/* wavenumber filter 			*/

	int ix1,ix2;		/* sample indices			*/
	int nx1fft,nx2fft;	/* dimensions after padding for FFT	*/
	int nK1,nK2;		/* transform (output) dimensions	*/
	float onK1,onK2;	/* 1 over transform (output) dimensions	*/
	float dK1,dK2;		/* wavenumber sampling interval		*/

        int ik1;		/* k1 counter				*/
        int ik2;		/* k2 counter				*/

	register complex **ct;	/* complex FFT workspace		*/
	register float **rt;	/* float FFT workspace			*/
	FILE *tracefp;		/* temp file to hold traces		*/
	FILE *hfp;		/* temp file to hold trace headers	*/

	float power;		/* power of i |k|			*/
	float arg;		/* argument of power			*/
	float sign;		/* sign of transform			*/
	float phasefac;		/* phase factor				*/
	float phase;		/* phase=phasefac*PI			*/
	complex exparg;		/* cwp_cexp(I arg)				*/

	cwp_Bool is_angle=cwp_false;	/* is angle set 		*/
	float angle=0.0;	/* angle for directional derivative	*/

	/* Hook up getpar to handle the parameters */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sukfrac );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get info from first trace */ 
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	if (tr.trid != TRID_DEPTH) warn("tr.trid = %d",tr.trid);
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


	if (!parObj.getparint("verbose",&verbose))	verbose = 1;
	if (!parObj.getparfloat("power",&power))	power = 0.0;
	if (!parObj.getparfloat("sign",&sign))		sign = 1.0;
	if (!parObj.getparfloat("phasefac",&phasefac))	phasefac = 0.0;
	phase = phasefac * PI;

	/* if angle is set then set is_angle boolean to true */
	if (parObj.getparfloat("angle",&angle)) is_angle = cwp_true;

	if (is_angle) {  /* convert angle to radians */
		angle = angle * PI/180.0;
	}
		

        /* Set up args for complex power evaluation */
        arg = sign * PIBY2 * power + phase;
        exparg = cwp_cexp(crmul(I, arg)); 

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
	if (nx2fft >= SU_NFLTS || nx2fft >= PFA_MAX)
		throw cseis_geolib::csException("Padded nx2=%d--too big",nx2fft);

	/* Determine number of wavenumbers in K1 and K2 */
	nK1 = nx1fft/2 + 1;
	nK2 = nx2fft/2 + 1;
	onK1 = 1.0/(float) nK1;
	onK2 = 1.0/(float) nK2;
	dK1 = TWOPI * onK1 / dx1;
	dK2 = TWOPI * onK2 / dx2;

	/* Allocate space */
	rt = ealloc2float(nx1fft, nx2fft);
	ct = ealloc2complex(nK1,nx2fft);
	kfilter = ealloc2complex(nx1fft,nx2fft);

	/* Zero all arrays */
	memset((void *) rt[0], 0, nx1fft*nx2fft*FSIZE);
	memset((void *) ct[0], 0, nK1*nx2fft*sizeof(complex));
	memset((void *) kfilter[0], 0, nx1fft*nx2fft*sizeof(complex));

	/* Build K-filter, by quadrant */
	/* positive k1, positive k2 */
	for (ik2=0; ik2<nK2; ++ik2) {
		for (ik1=0; ik1<nK1; ++ik1)  {

			/* find the magnitude of the k-vector */
			/* Build filter */
			if (is_angle) {
				k = sqrt(ik1*ik1*dK1*dK1*cos(angle)*cos(angle) + ik2*ik2*dK2*dK2*sin(angle)*sin(angle));
			} else {

				k = sqrt(ik1*ik1*dK1*dK1 + ik2*ik2*dK2*dK2);
			}

			/* kludge to handle k=0 case */
			if (power < 0 && ik1==0 && ik2==0 ) k = FLT_MAX;

			/* Calculate filter value */
			amp = pow(k,power)* onK1 * onK2 /PI;
			kfilter[ik2][ik1] = crmul(exparg, amp);
		}
	}


	/* positive k1, negative k2 */
	for (ik2=nK2; ik2<nx2fft; ++ik2) {
		for (ik1=0 ; ik1<nK1; ++ik1) {

			/* find the magnitude of the k-vector */
			if (is_angle) {
				k = sqrt(ik1*ik1*dK1*dK1*cos(angle)*cos(angle) +
					  (nx2fft-ik2)*(nx2fft -ik2)*dK2*dK2*sin(angle)*sin(angle) );

			} else {

				k = sqrt(ik1*ik1*dK1*dK1 +
					  (nx2fft-ik2)*(nx2fft -ik2)*dK2*dK2 );
			}
			/* Calculate filter value */
			amp = pow(k,power) * onK1 * onK2/PI;

			/* Build filter */
			kfilter[ik2][ik1] = crmul(exparg, amp);
		}
	}

	/* Load traces into fft arrays and close tmpfile */
	erewind(tracefp);
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
			ct[ik2][ik1] = cmul(ct[ik2][ik1], kfilter[ik2][ik1]);


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


} // END namespace
