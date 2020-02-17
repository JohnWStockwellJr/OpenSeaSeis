/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SULFAF: $Revision: 1.3 $ ; $Date: 2011/11/12 00:09:00 $	*/

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
std::string sdoc_sulfaf =
" SULFAF -  Low frequency array forming					" 
"			  						"
"  sulfaf < stdin > stdout [optional parameters]			"
"			  						"
" Optional Parameters:	  						"
" key=ep	header keyword describing the gathers			"
" f1=3		lower frequency	cutof					"
" f2=20		high frequency cutof					"
" fr=5		frequency ramp						"
"			  						"
" vel=1000	surface wave velocity					"
" dx=10		trace spacing						"
" maxmix=tr.ntr	default is the entire gather				"
" adb=1.0	add back ratio 1.0=pure filtered 0.0=origibal		"
"			  						"
" tri=0		1 Triangular weight in mixing window			"
"			  						"
" Notes:		  						"
" The traces transformed into the freqiency domain			"
" where a trace mix is performed in the specifed frequency range	"
" as Mix=vel/(freq*dx)							"

" This program uses \"get_gather\" and \"put_gather\" so requires that  "
" the  data be sorted into ensembles designated by \"key\" with the ntr"
" field set to the number of traces in each respective ensemble.	"
"									"
" Example:								 "
" susort ep offset < data.su > datasorted.su				"
" suputgthr dir=Data verbose=1 < datasorted.su			  "
" sugetgthr dir=Data verbose=1 > dataupdated.su			 "
" sulfaf  < dataupdated.su > ccfiltdata.su				"
"			  						"
" (Work in progress, editing required)                 			"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sulfaf {


#define LOOKFAC 1	/* Look ahead factor for npfaro */
#define PFA_MAX 720720  /* Largest allowed nfft		*/
#define PIP2 PI/2.0	/* IP/2				*/

void* main_sulfaf( void* args )
{
	cwp_String key;		/* header key word from segy.h		*/
	cwp_String type;	/* ... its type				*/
	Value val;		/* ... its value			*/
	segy **rec_o;		/* trace header+data matrix		*/
	int first=0;		/* true when we passed the first gather */
	int ng=0;
	float dt;		/* time sampling interval	*/
	int nt;			/* num time samples per trace	*/
	int ntr;		/* num of traces per ensemble	*/
	
	int nfft=0;		/* lenghth of padded array	*/
	float snfft;		/* scale factor for inverse fft */
	int nf=0;		/* number of frequencies	*/
	float d1;		/* frequency sampling interval.	*/
	float *rt=NULL;		/* real trace			*/
	complex *ct=NULL;	/* complex trace		*/
	float **ffdr=NULL;	/* frequency domain data	*/
	float **ffdi=NULL;	/* frequency domain data	*/
	float **ffdrm=NULL;	/* frequency domain mixed data	*/
	float **ffdim=NULL;	/* frequency domain mixed data	*/
	
	int verbose;		/* flag: =0 silent; =1 chatty	*/
	
	
	float f1;		/* minimum frequency		*/
	int if1;		/* ...    ...  integerized	*/
	float f2;		/* maximum frequency		*/
	int if2;		/* ...    ...  integerized	*/
	float fr;		/* slope of frequency ramp	*/
	int ifr;		/* ...    ...  integerized	*/
	float vel;		/* velocity of guided waves	*/
	float dx;		/* spatial sampling intervall	*/
	int maxmix;		/* size of mix			*/
	int tri;		/* flag: =1 trianglular window	*/
	float adb;		/* add back ratio		*/
		
	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sulfaf );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	
	if (!parObj.getparstring("key", &key))		key = "ep";
	if (!parObj.getparfloat("f1", &f1))		f1 = 3.0;
	if (!parObj.getparfloat("f2", &f2))		f2 = 20.0;
	if (!parObj.getparfloat("dx", &dx))		dx = 10;
	if (!parObj.getparfloat("vel", &vel))		vel = 1000;
	if (!parObj.getparfloat("fr", &fr))		fr = 5;
	if (!parObj.getparint("maxmix", &maxmix))	maxmix = -1;
	if (!parObj.getparint("tri", &tri))		tri = 0;
	if (!parObj.getparfloat("adb", &adb))		adb = 1.0;
	
	if (!parObj.getparint("verbose", &verbose)) verbose = 0;

	/* get the first record */
	rec_o = get_gather(&key,&type,&val,&nt,&ntr,&dt,&first);
	if(ntr==0) throw cseis_geolib::csException("Can't get first record\n");
		
	/* set up the fft */
	nfft = npfaro(nt, LOOKFAC * nt);
	 if (nfft >= SU_NFLTS || nfft >= PFA_MAX)
		 	throw cseis_geolib::csException("Padded nt=%d--too big", nfft);
	 nf = nfft/2 + 1;
	 snfft=1.0/nfft;
	d1=1.0/(nfft*dt);
	
	ct=ealloc1complex(nf);
	rt=ealloc1float(nfft);
	
	if1=NINT(f1/d1);
	if2=NINT(f2/d1);
	ifr=NINT(fr/d1);
	
	do {
		if(maxmix==-1) maxmix=ntr;
		ng++;
		
		/* Allocate arrays for fft */
		ffdr = ealloc2float(nf,ntr);
		ffdi = ealloc2float(nf,ntr);
		ffdrm = ealloc2float(if2+ifr,ntr);
		ffdim = ealloc2float(if2+ifr,ntr);
		{ int itr,iw;
			for(itr=0;itr<ntr;itr++) {
				
				memcpy( (void *) rt, (const void *) (*rec_o[itr]).data, nt*FSIZE);

				memset( (void *) &rt[nt], 0,(nfft-nt)*FSIZE);
				
				pfarc(1,nfft,rt,ct);
				
				for(iw=0;iw<nf;iw++) {
					ffdr[itr][iw] = ct[iw].r;
					ffdi[itr][iw] = ct[iw].i;
				}
			}
		}
		
		/* Mixing */
		{ int mix,iw,nmix;
		  int ims,ime;
		  int itr,itrm,iww,ws;
		  float tmpr,tmpi,wh;
		
			for(iw=if1;iw<if2+ifr;iw++) {
				
				mix=MIN(NINT(vel/iw*d1*dx),maxmix);
				if(!ISODD(mix)) mix -=1;
				if (verbose) warn(" %f %d",iw*d1,mix);
				
				for(itr=0;itr<ntr;itr++) {
						
					ims=MAX(itr-mix/2,0);
					ime=MIN(ims+mix,ntr-1);

					tmpr=0.0; tmpi=0.0;
					wh=1.0; ws=mix/2;
					nmix=0;
					for(itrm=ims,iww=-mix/2;itrm<ime;++itrm,++iww) {
						++nmix;
						if(tri) wh = (float)(ws-abs(iww));
						tmpr+=ffdr[itrm][iw]*wh;
						tmpi+=ffdi[itrm][iw]*wh;
					}
					ffdrm[itr][iw]=tmpr/nmix;
					ffdim[itr][iw]=tmpi/nmix;
				}
			}
			
			for(iw=if1;iw<if2;iw++) {
				for(itr=0;itr<ntr;itr++) {
					ffdr[itr][iw]=ffdrm[itr][iw]*adb+ffdr[itr][iw]*(1.0-adb);
					ffdi[itr][iw]=ffdim[itr][iw]*adb+ffdi[itr][iw]*(1.0-adb);
				}
			}
			
			for(iw=if2,iww=0;iw<if2+ifr;iw++,iww++) {
				
				wh=(float)(1.0-(float)iww/(float)ifr);
				
				for(itr=0;itr<ntr;itr++) {
					ffdr[itr][iw] = (wh*ffdrm[itr][iw]+ffdr[itr][iw]*(1.0-wh))*adb+ffdr[itr][iw]*(1.0-adb);
					ffdi[itr][iw] = (wh*ffdim[itr][iw]+ffdi[itr][iw]*(1.0-wh))*adb+ffdi[itr][iw]*(1.0-adb);
				}
			}
		
		}
		  
		
		{ int itr,iw;
			for(itr=0;itr<ntr;itr++) {
				
				for(iw=0;iw<nf;iw++) {
					ct[iw].r = ffdr[itr][iw]*snfft;
					ct[iw].i = ffdi[itr][iw]*snfft;
				} 
				
				pfacr(-1,nfft,ct,rt);
				memcpy( (void *) (*rec_o[itr]).data, (const void *) rt, nt*FSIZE);
				
			}
		}
		
			rec_o = put_gather(rec_o,&nt,&ntr);

		free2float(ffdr);
		free2float(ffdi);
		free2float(ffdrm);
		free2float(ffdim);
		rec_o = get_gather(&key,&type,&val,&nt,&ntr,&dt,&first);
		
	} while(ntr);

	warn("Number of gathers %10d\n",ng);
	 
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
