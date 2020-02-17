/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSMGAUSS2: $Revision: 1.2 $ ; $Date: 2011/11/12 00:09:00 $	*/

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
std::string sdoc_susmgauss2 =
"									"
" SUSMGAUSS2 --- SMOOTH a uniformly sampled 2d array of velocities	"
"		using a Gaussian filter specified with correlation 	"
"	lengths a1 and a2.	"
"									"
" susmgauss2 < stdin [optional parameters ] > stdout			"
"									"
" Optional Parameters:							"
" a1=0			smoothing parameter in the 1 direction		"
" a2=0			smoothing parameter in the 2 direction		"
"									"
" Notes:								"
" Larger a1 and a2 result in a smoother velocity. The velocities are	"
" first transformed to slowness and then a Gaussian filter is applied	"
" in the wavenumber domain.						"
"	"
" Input file must be in SU format. The output file is smoothed velocity "

;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace susmgauss2 {


/*
 *	Credits: 
 *		CWP: Carlos Pacheco, 2005
 */
/**************** end self doc *******************************************/

/* definitions */
#define PFA_MAX 720720   /* Largest allowed nfft */
#define LOOKFAC 2	/* Look factor */

segy tr;

void* main_susmgauss2( void* args )
{
	int nx1,nx2;		/* numbers of samples			*/
	int ix1, ix2;		/* sample indices			*/
	float a1,a2;		/* filter dimensions			*/

	float nyq;		/* Nyquist wavenumber			*/
	float pi;		/* pi number				*/
	float vmax;		/* maximum value of the data		*/
	float vfmax;		/* scale factor after filtering		*/
	float c1,c2;
	float **v=NULL;		/* array of velocities			*/
	float *k1=NULL,*k2=NULL;	/* wavenumber arrays		*/
	float *kfilt1=NULL,*kfilt2=NULL;/* intermediate filter arrays   */
	float dk1,dk2;			/* wavenumber interval		*/
	float **kfilter=NULL;		/* array of filter values	*/

	int nx1fft,nx2fft;	/* dimensions after padding for FFT	*/
	int nK1,nK2;	  /* transform dimension			*/
	int ik1,ik2;	  /* wavenumber indices				*/

	register complex **ct=NULL;	/* complex FFT workspace	*/
	register float **rt=NULL;	/* float FFT workspace		*/
	FILE *tracefp=NULL;	/* temp file to hold traces		*/
	FILE *hfp=NULL;		/* temp file to hold trace headers	*/
	
	/* hook up getpar to handle the parameters */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_susmgauss2 );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get parameters from command line */
	if (!parObj.getparfloat("a1",&a1)) a1=0.;
	if (!parObj.getparfloat("a2",&a2)) a2=0.;

	/* Get info from first trace */
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	if (tr.trid != TRID_DEPTH)  warn("tr.trid=%d",tr.trid);
	nx1=tr.ns;

	/* Compute Nyquist wavenumber */
	nyq=0.5;
	
	/* Store traces in tmpfile while getting a count */
	tracefp=etmpfile();
	hfp=etmpfile();
	nx2=0;
	do {
		++nx2;
		efwrite(&tr,HDRBYTES, 1, hfp);
		efwrite(tr.data, FSIZE, nx1, tracefp);
	} while (cs2su->getTrace(&tr));
	
	/* Determine number of wavenumbers in K1 and K2 */
	nx1fft=npfaro(nx1, LOOKFAC*nx1);
	nx2fft=npfa(nx2);
	if (nx1fft >=SU_NFLTS || nx1fft >= PFA_MAX)
	  throw cseis_geolib::csException("Padded nx1=%d--too big",nx1fft);
	if (nx2fft >= PFA_MAX)
	  throw cseis_geolib::csException("Padded nx2=%d--too big",nx2fft);

	/* Determine number of wavenumbers in K1 and K2 */
	nK1=nx1fft/2 + 1;
	nK2=nx2fft/2 + 1;

	/* Allocate space */
	v=alloc2float(nx1,nx2);
	rt=alloc2float(nx1fft,nx2fft);
	ct=alloc2complex(nK1,nx2fft);
	kfilter=alloc2float(nx1fft,nx2fft);
	k1=alloc1float(nK1);
	k2=alloc1float(nK2);
	kfilt1= alloc1float(nK1);
	kfilt2= alloc1float(nK2);

	/* Zero all arrays */
	memset((void *) rt[0], 0, nx1fft*nx2fft*FSIZE);
	memset((void *) kfilter[0], 0, nx1fft*nx2fft*FSIZE);
	memset((void *) ct[0], 0, nK1*nx2fft*sizeof(complex));
	memset((void *) k1, 0, nK1*FSIZE);
	memset((void *) k2, 0, nK2*FSIZE);
	memset((void *) kfilt1, 0,  nK1*FSIZE);
	memset((void *) kfilt2, 0, nK2*FSIZE);

	/* Determine wavenumber arrays for the filter */
	pi=PI;
	dk1=2*pi / nx1fft;
	dk2=2*pi / nx2fft;
	for (ik1=0; ik1<nK1; ++ik1) {
	  c1=a1*ik1*dk1/ 2;
	  kfilt1[ik1]= exp(-pow(c1,2.0));
	}
	for (ik2=0; ik2<nK2; ++ik2) {
	  c2= a2*ik2*dk2/2;
	  kfilt2[ik2]= exp(-pow(c2,2.0));
	}
	
	/* Build Gaussian filter */
	/* positive k1, positive k2 */
	for (ik2=0; ik2<nK2; ++ik2) {
	  for (ik1=0; ik1<nK1; ++ik1) {
		kfilter[ik2][ik1]=kfilt2[ik2]*kfilt1[ik1];
	  }
	}
	/* positive k1, negative k2 */
	for (ik2=nK2; ik2<nx2fft; ++ik2) {
	  for (ik1=0; ik1<nK1; ++ik1) {
		kfilter[ik2][ik1]=kfilt2[nx2fft-ik2]*kfilt1[ik1];
	  }
	}

	/* Read velocities from temp file and determine maximum */
	rewind(tracefp);
	fread(v[0],sizeof(float),nx2*nx1,tracefp);
	vmax=v[0][0];
	for (ix2=0; ix2<nx2; ++ix2) {
	  for (ix1=0; ix1<nx1; ++ix1) {
		vmax=MAX(vmax,v[ix2][ix1]);
	  }
	}

	/* Load data into FFT arrays */
	rewind(tracefp);
	for (ix2=0; ix2<nx2; ++ix2) {
	  efread(rt[ix2], FSIZE, nx1, tracefp);
	}   
	
	/* Fourier transform dimension 1 */
	pfa2rc(-1,1,nx1fft,nx2,rt[0],ct[0]);
	
	/* Fourier transform dimension 2 */
	pfa2cc(-1,2,nK1,nx2fft,ct[0]);

	/* Apply filter to the data */
	for (ik2=0; ik2<nx2fft; ++ik2) {
	  	for (ik1=0; ik1<nK1;  ++ik1) {
			ct[ik2][ik1]=crmul(ct[ik2][ik1], kfilter[ik2][ik1]) ;
		}
	}

	/* Inverse Fourier transformation dimension 2 */
	pfa2cc(1,2,nK1,nx2fft,ct[0]);

	/* Inverse Fourier transformation dimension 1 */
	pfa2cr(1,1,nx1fft,nx2,ct[0],rt[0]);

	/* Find maximum of filtered data */
	vfmax=rt[0][0];
	for (ix2=0; ix2<nx2; ++ix2) {
		for (ix1=0; ix1<nx1; ++ix1) {
			vfmax=MAX(vfmax,rt[ix2][ix1]);
		}
	}

	/* Rescale and output filtered data */
	erewind(hfp);
	for (ix2=0; ix2<nx2; ++ix2) {
		efread(&tr, HDRBYTES, 1, hfp);
		for (ix1=0; ix1<nx1; ++ix1)
			tr.data[ix1]=(rt[ix2][ix1]) * vmax / vfmax;
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
