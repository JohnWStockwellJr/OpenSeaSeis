/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.		       */

/* SUKTMIG2D: $Revision: 1.6 $ ; $Date: 2011/12/23 19:33:01 $*/

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
std::string sdoc_suktmig2d =
" 								"
" SUKTMIG2D - prestack time migration of a common-offset	"
"	section with the double-square root (DSR) operator	"
"								"
" 								"
"   suktmig2d < infile vfile= [parameters]  > outfile		"
" 								"
" Required Parameters:						"
" vfile=	rms velocity file (units/s) v(t,x) as a function"
"		of time						"
" dx=		distance (units) between consecutive traces	"
"								"
" Optional parameters:						"
" fcdpdata=tr.cdp	first cdp in data			"
" firstcdp=fcdpdata	first cdp number in velocity file	"
" lastcdp=from header	last cdp number in velocity file	"
" dcdp=from header	number of cdps between consecutive traces"
" angmax=40	maximum aperture angle for migration (degrees)	"
" hoffset=.5*tr.offset		half offset (m)			"
" nfc=16	number of Fourier-coefficients to approximate	"
"		low-pass					"
"		filters. The larger nfc the narrower the filter	"
" fwidth=5 	high-end frequency increment for the low-pass	"
" 		filters						"
" 		in Hz. The lower this number the more the number"
"		of lowpass filters to be calculated for each 	"
"		input trace.					"
"								"
" Caveat: this code may need some work				"
" Notes:							"
" Data must be preprocessed with sufrac to correct for the	"
" wave-shaping factor using phasefac=.25 for 2D migration.	"
"								"
" Input traces must be sorted into offset and cdp number.	"
" The velocity file consists of rms velocities for all CMPs as a"
" function of vertical time and horizontal position v(t,z)	"
" in C-style binary floating point numbers. It's easiest to 	"
" supply v(t,z) that has the same dimensions as the input data to"
" be migrated.							"
" The units may be feet or meters, as long as these are		"
" consistent.							"
" Antialias filter is performed using (Gray,1992, Geoph. Prosp), "
" using nc low- pass filtered copies of the data. The cutoff	"
" frequencies are calculated  as fractions of the Nyquist	"
" frequency.							"
"								"
" The maximum allowed angle is 80 degrees(a 10 degree taper is "
" applied to the end of the aperture)				"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suktmig2d {


#define LOOKFAC 2       /* Look ahead factor for npfaro   */
#define PFA_MAX 720720  /* Largest allowed nfft	   */


/* Prototype of functions used internally */
void lpfilt(int nfc, int nfft, float dt, float fhi, float *filter);

segy intrace; 	/* input traces */
segy outtrace;	/* migrated output traces */

void* main_suktmig2d( void* args )
{
	int i,k,imp,iip,it,ix,ifc;	/* counters */
	int ntr,nt;			/* x,t */

	int verbose;	/* is verbose?				*/
	int nc;		/* number of low-pass filtered versions	*/
			/*  of the data for antialiasing	*/
	int nfft,nf;	/* number of frequencies		*/
	int nfc;	/* number of Fourier coefficients for low-pass filter */
	int fwidth;	/* high-end frequency increment for the low-pass */
			/* filters 				*/
	int firstcdp=0;	/* first cdp in velocity file		*/
	int lastcdp=0;	/* last cdp in velocity file		*/
	int oldcdp=0;	/* temporary storage			*/
	int fcdpdata=0;	/* first cdp in the data		*/
	int olddeltacdp=0;
	int deltacdp;
	int ncdp=0;	/* number of cdps in the velocity file	*/
	int dcdp=0;	/* number of cdps between consecutive traces */

	float dx=0.0;	/* cdp sample interval */
	float hoffset=0.0;  /* half receiver-source */
	float p=0.0;	/* horizontal slowness of the migration operator */
	float pmin=0.0;	/* maximum horizontal slowness for which there's */
			/* no aliasing of the operator */
	float dt;	/* t sample interval */
	float h;	/* offset */
	float x;	/* aperture distance */
	float xmax=0.0;	/* maximum aperture distance */

	float obliq;	/* obliquity factor */
	float geoms;	/* geometrical spreading factor */
	float angmax;   /* maximum aperture angle */

	float mp,ip;	/* mid-point and image-point coordinates */
	float t;	/* time */
	float t0;	/* vertical traveltime */
	float tmax;	/* maximum time */

	float fnyq;	/* Nyquist frequency */
	float ang;	/* aperture angle */
	float angtaper=0.0;	/* aperture-angle taper */
	float v;		/* velocity */
  
	float *fc=NULL;		/* cut-frequencies for low-pass filters */
	float *filter=NULL;	/* array of low-pass filter values */

	float **vel=NULL;	/* array of velocity values from vfile */
	float **data=NULL;	/* input data array*/
	float **lowpass=NULL;   /* low-pass filtered version of the trace */
	float **mig=NULL;	/* output migrated data array */

	register float *rtin=NULL,*rtout=NULL;/* real traces */
	register complex *ct=NULL;   /* complex trace */

	/* file names */
	char *vfile="";		/* name of velocity file */
	FILE *vfp=NULL;
	FILE *tracefp=NULL;	/* temp file to hold traces*/
	FILE *hfp=NULL;		/* temp file to hold trace headers */

	float datalo[8], datahi[8];
	int itb, ite;
	float firstt, amplo, amphi;

	cwp_Bool check_cdp=cwp_false;	/* check cdp in velocity file	*/

	/* Hook up getpar to handle the parameters */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suktmig2d );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	
	/* Get info from first trace */
	if (!cs2su->getTrace(&intrace))  throw cseis_geolib::csException("can't get first trace");
	nt=intrace.ns;
	dt=(float)intrace.dt/1000000;
	tmax=(nt-1)*dt;

	CSMUSTGETPARFLOAT("dx",&dx);
	CSMUSTGETPARSTRING("vfile",&vfile);
	if (!parObj.getparfloat("angmax",&angmax)) angmax=40;
	if (!parObj.getparint("firstcdp",&firstcdp)) firstcdp=intrace.cdp;
	if (!parObj.getparint("fcdpdata",&fcdpdata)) fcdpdata=intrace.cdp;
	if (!parObj.getparfloat("hoffset",&hoffset)) hoffset=.5*intrace.offset;
	if (!parObj.getparint("nfc",&nfc)) nfc=16;
	if (!parObj.getparint("fwidth",&fwidth)) fwidth=5;
	if (!parObj.getparint("verbose",&verbose)) verbose=0;

	h=hoffset;

	/* Store traces in tmpfile while getting a count of number of traces */
	tracefp = etmpfile();
	hfp = etmpfile();
	ntr = 0;
	do {
		++ntr;

		/* get new deltacdp value */
		deltacdp=intrace.cdp-oldcdp;

		/* read headers and data */
		efwrite(&intrace,HDRBYTES, 1, hfp);
		efwrite(intrace.data, FSIZE, nt, tracefp);

		/* error trappings. */
		/* ...did cdp value interval change? */
		if ((ntr>3) && (olddeltacdp!=deltacdp)) {

			if (verbose) {
			warn("cdp interval changed in data");	
			warn("ntr=%d olddeltacdp=%d deltacdp=%d"
				,ntr,olddeltacdp,deltacdp);
		 	check_cdp=cwp_true;
			}
		}
		
		/* save cdp and deltacdp values */
		oldcdp=intrace.cdp;
		olddeltacdp=deltacdp;

	} while (cs2su->getTrace(&intrace));

	/* get last cdp  and dcdp */
	if (!parObj.getparint("lastcdp",&lastcdp)) lastcdp=intrace.cdp; 
	if (!parObj.getparint("dcdp",&dcdp))	dcdp=deltacdp - 1;


	parObj.checkpars();

	/* error trappings */
	if ( (firstcdp==lastcdp) 
		|| (dcdp==0) 
		|| (check_cdp==cwp_true) ) warn("Check cdp values in data!");

	/* rewind trace file pointer and header file pointer */
	erewind(tracefp);
	erewind(hfp);

	/* total number of cdp's in data */
	ncdp=lastcdp-firstcdp+1;

	/* Set up FFT parameters */
	nfft = npfaro(nt, LOOKFAC*nt);
	if(nfft>= SU_NFLTS || nfft >= PFA_MAX)
	  throw cseis_geolib::csException("Padded nt=%d -- too big",nfft);
	nf = nfft/2 + 1;

	/* Determine number of filters for antialiasing */
	fnyq= 1.0/(2*dt);
	nc=ceil(fnyq/fwidth);
	if (verbose)
		warn(" The number of filters for antialiasing is nc= %d",nc);

	/* Allocate space */
	data = alloc2float(nt,ntr);
	lowpass=alloc2float(nt,nc+1);
	mig=   alloc2float(nt,ntr);
	vel=   alloc2float(nt,ncdp);
	fc = alloc1float(nc+1);
	rtin= ealloc1float(nfft);
	rtout= ealloc1float(nfft);
	ct= ealloc1complex(nf);
	filter= alloc1float(nf);

	/* Read data from temporal array */
	for (ix=0; ix<ntr; ++ix){
		efread(data[ix],FSIZE,nt,tracefp);
	}

	/* read velocities */
	vfp=efopen(vfile,"r");
	efread(vel[0],FSIZE,nt*ncdp,vfp);
	efclose(vfp);

	/* Zero all arrays */
	memset((void *) mig[0], 0,nt*ntr*FSIZE);
	memset((void *) rtin, 0, nfft*FSIZE);
	memset((void *) filter, 0, nf*FSIZE);
	memset((void *) lowpass[0], 0,nt*(nc+1)*FSIZE);

	/* Calculate cut frequencies for low-pass filters */
	for(i=1; i<nc+1; ++i){
		fc[i]= fnyq*i/nc;
	}

	/* Start the migration process */
	/* Loop over input mid-points first */
	if (verbose) warn("Starting migration process...\n");
	for(imp=0; imp<ntr; ++imp){
		float perc;

		mp=imp*dx; 
		perc=imp*100.0/(ntr-1);
		if(fmod(imp*100,ntr-1)==0 && verbose)
			warn("migrated %g\n ",perc);

		/* Calculate low-pass filtered versions  */
		/* of the data to be used for antialiasing */
		for(it=0; it<nt; ++it){
			rtin[it]=data[imp][it];
		}
		for(ifc=1; ifc<nc+1; ++ifc){
			memset((void *) rtout, 0, nfft*FSIZE);
			memset((void *) ct, 0, nf*FSIZE);
			lpfilt(nfc,nfft,dt,fc[ifc],filter);
			pfarc(1,nfft,rtin,ct);

			for(it=0; it<nf; ++it){
				ct[it] = crmul(ct[it],filter[it]);
			}
			pfacr(-1,nfft,ct,rtout);
			for(it=0; it<nt; ++it){ 
				lowpass[ifc][it]= rtout[it]; 
			}
		}

		/* Loop over vertical traveltimes */
		for(it=0; it<nt; ++it){
			int lx,ux;

			t0=it*dt;
			v=vel[imp*dcdp+fcdpdata-1][it];
			xmax=tan((angmax+10.0)*PI/180.0)*v*t0;
			lx=MAX(0,imp - ceil(xmax/dx)); 
			ux=MIN(ntr,imp + ceil(xmax/dx));
	
		/* loop over output image-points to the left of the midpoint */
		for(iip=imp; iip>lx; --iip){
			float ts,tr;
			int fplo=0, fphi=0;
			float ref,wlo,whi;

			ip=iip*dx; 
			x=ip-mp; 
			ts=sqrt( pow(t0/2,2.0) + pow((x+h)/v,2.0) );
			tr=sqrt( pow(t0/2,2.0) + pow((h-x)/v,2.0) );
			t= ts + tr;
			if(t>=tmax) break;
			geoms=sqrt(1/(t*v));
	  		obliq=sqrt(.5*(1 + (t0*t0/(4*ts*tr)) 
					- (1/(ts*tr))*sqrt(ts*ts - t0*t0/4)*sqrt(tr*tr - t0*t0/4)));
	  		ang=180.0*fabs(acos(t0/t))/PI;  
	  		if(ang<=angmax) angtaper=1.0;
	  		if(ang>angmax) angtaper=cos((ang-angmax)*PI/20);
	  		/* Evaluate migration operator slowness p to determine */
			/* the low-pass filtered trace for antialiasing */
			pmin=1/(2*dx*fnyq);
			p=fabs((x+h)/(pow(v,2.0)*ts) + (x-h)/(pow(v,2.0)*tr));
				if(p>0){fplo=floor(nc*pmin/p);}
				if(p==0){fplo=nc;}
				ref=fmod(nc*pmin,p);
				wlo=1-ref;
				fphi=++fplo;
				whi=ref;
				itb=MAX(ceil(t/dt)-3,0);
				ite=MIN(itb+8,nt);
				firstt=(itb-1)*dt;
				/* Move energy from CMP to CIP */
				if(fplo>=nc){
					for(k=itb; k<ite; ++k){
						datalo[k-itb]=lowpass[nc][k];
					}
					ints8r(8,dt,firstt,datalo,0.0,0.0,1,&t,&amplo);
					mig[iip][it] +=geoms*obliq*angtaper*amplo;
				} else if(fplo<nc){
					for(k=itb; k<ite; ++k){
						datalo[k-itb]=lowpass[fplo][k];
						datahi[k-itb]=lowpass[fphi][k];
					}
					ints8r(8,dt,firstt,datalo,0.0,0.0,1,&t,&amplo);
					ints8r(8,dt,firstt,datahi,0.0,0.0,1,&t,&amphi);
					mig[iip][it] += geoms*obliq*angtaper*(wlo*amplo + whi*amphi);
				}
			}

			/* loop over output image-points to the right of the midpoint */
			for(iip=imp+1; iip<ux; ++iip){
				float ts,tr;
				int fplo=0, fphi;
				float ref,wlo,whi;

				ip=iip*dx; 
				x=ip-mp; 
				t0=it*dt;	  
				ts=sqrt( pow(t0/2,2.0) + pow((x+h)/v,2.0) );
				tr=sqrt( pow(t0/2,2.0) + pow((h-x)/v,2.0) );
				t= ts + tr;
				if(t>=tmax) break;
				geoms=sqrt(1/(t*v));
				obliq=sqrt(.5*(1 + (t0*t0/(4*ts*tr)) 
					- (1/(ts*tr))*sqrt(ts*ts 
						- t0*t0/4)*sqrt(tr*tr 
								- t0*t0/4)));
				ang=180.0*fabs(acos(t0/t))/PI;   
				if(ang<=angmax) angtaper=1.0;
				if(ang>angmax) angtaper=cos((ang-angmax)*PI/20.0);

				/* Evaluate migration operator slowness p to determine the  */
				/* low-pass filtered trace for antialiasing */
				pmin=1/(2*dx*fnyq);
				p=fabs((x+h)/(pow(v,2.0)*ts) + (x-h)/(pow(v,2.0)*tr));
				if(p>0){
					fplo=floor(nc*pmin/p);
				}
				if(p==0){
					fplo=nc;
				}

				ref=fmod(nc*pmin,p);
				wlo=1-ref;
				fphi=fplo+1;
				whi=ref;
				itb=MAX(ceil(t/dt)-3,0);
				ite=MIN(itb+8,nt);
				firstt=(itb-1)*dt;

				/* Move energy from CMP to CIP */
				if(fplo>=nc){
					for(k=itb; k<ite; ++k){
						datalo[k-itb]=lowpass[nc][k];
					}
					ints8r(8,dt,firstt,datalo,0.0,0.0,1,&t,&amplo);
					mig[iip][it] +=geoms*obliq*angtaper*amplo;
				} else if(fplo<nc){
					for(k=itb; k<ite; ++k){
						datalo[k-itb]=lowpass[fplo][k];
						datahi[k-itb]=lowpass[fphi][k];
					}
					ints8r(8,dt,firstt,datalo,0.0,0.0,1,&t,&amplo);
					ints8r(8,dt,firstt,datahi,0.0,0.0,1,&t,&amphi);
					mig[iip][it] += geoms*obliq*angtaper*(wlo*amplo + whi*amphi);
				}
			}

		}
	} 

	/* Output migrated data */
	erewind(hfp);
	for (ix=0; ix<ntr; ++ix) {
		efread(&outtrace, HDRBYTES, 1, hfp);
		for (it=0; it<nt; ++it) {
			outtrace.data[it] = mig[ix][it];
		}
		su2cs->putTrace(&outtrace);
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

void
lpfilt(int nfc, int nfft, float dt, float fhi, float *filter)
/*******************************************************************************
lpfilt -- low-pass filter using Lanczos Smoothing 
	(R.W. Hamming:"Digital Filtering",1977)
****************************************************************************
Input: 
nfc	number of Fourier coefficients to approximate ideal filter
nfft	number of points in the fft
dt	time sampling interval
fhi	cut-frequency

Output:
filter  array[nf] of filter values
*****************************************************************************
Notes: Filter is to be applied in the frequency domain   
*****************************************************************************
Author: CWP: Carlos Pacheco   2006   
*****************************************************************************/
{
	int i,j;  /* counters */
	int nf;   /* Number of frequencies (including Nyquist) */
	float onfft;  /* reciprocal of nfft */
	float fn; /* Nyquist frequency */
	float df; /* frequency interval */
	float dw; /* frequency interval in radians */
	float whi;/* cut-frequency in radians */
	float w;  /* radian frequency */

	nf= nfft/2 + 1;
	onfft=1.0/nfft;
	fn=1.0/(2*dt);
	df=onfft/dt;
	whi=fhi*PI/fn;
	dw=df*PI/fn;

	for(i=0; i<nf; ++i){
		filter[i]= whi/PI;
		w=i*dw;

		for(j=1; j<nfc; ++j){
			float c= sin(whi*j)*sin(PI*j/nfc)*2*nfc/(PI*PI*j*j);
			filter[i] +=c*cos(j*w);
		}
	}
}
  

  
  
  

  
  
  

} // END namespace
