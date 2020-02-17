/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUFWATRIM: $Revision: 1.2 $ ; $Date: 2011/11/12 00:09:00 $	  */


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
std::string sdoc_sufwatrim =
" SUFWATRIM - FX domain Alpha TRIM					"
"			 						"
"  sufwatrim  <stdin > stdout [optional parameters]			"
"			 						"
" Required parameters:							"
" key=key1,key2,..	Header words defining mixing dimesnion		"
" dx=d1,d2,..		Distance units for each header word		"
" Optional parameters:							"
" keyg=ep		Header word indicating the start of gather	"
" vf=0			=1 Do a frequency dependent mix			"
" vmin=5000		Velocity of the reflection slope		"
"			than should not be attenuated			"
" Notes:		 						"
" Trace with the header word mark set to one will be 			"
" the output trace 							"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sufwatrim {


/*
 * Credits: Potash Corporation of Saskatchewan, Balasz Nemeth
 * Code given to CWP in 2008
 */


/**************** end self doc ********************************/

#define PFA_MAX 720720  /* Largest allowed nfft	   */
#define PIP2 1.570796327
#define AMPSP(c) rcabs(c)
#define PHSSP(c) atan2(c.i,c.r)

segy tr;
float n_distance(segy **rec_o,int *index,cwp_String *type,float *dx,unsigned int nd,unsigned int imx,unsigned int itr);
float alpha_trim(float *a,int n,float p);
float alpha_trim_w(float *a,float *w,int n,float p);

int verbose;

void* main_sufwatrim( void* args )
{
	cwp_String keyg;	/* header key word from segy.h		*/
	cwp_String typeg;	/* ... its type				*/
	Value valg;
	cwp_String key[SU_NKEYS];	/* array of keywords		 */
	cwp_String type[SU_NKEYS];	/* array of keywords		 */
	int index[SU_NKEYS];	/* name of type of getparred key	 */
	
	segy **rec_o;		/* trace header+data matrix */	
	
	int first=0;	/* true when we passed the first gather */
	int ng=0;
	float dt;	/* time sampling interval		*/
	int nt;		/* number of time samples per trace	*/
	int ntr;	/* number of traces per ensemble	*/
	
	int nfft=0;		/* lenghth of padded array		*/
	float snfft;		/* scale factor for inverse fft		*/
	int nf=0;		/* number of frequencies		*/
	float d1;		/* frequency sampling int.		*/
	float *rt;		/* real trace				*/
	complex *ctmix;		/* complex trace			*/
	complex **fd;		/* frequency domain data		*/

	
	float padd;
	
	int nd;			/* number of dimensions */
	float *dx=NULL;
	float fac;
	float vmin;
	int vf;
	
	/* Trimming arrays */
	float *itrm=NULL;
	float *rtrm=NULL;
	float *wht=NULL;
	float trimp=15;
		
	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sufwatrim );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	
	if (!parObj.getparstring("keyg", &keyg)) keyg ="ep";
	if (!parObj.getparint("vf", &vf)) vf = 1;
	if (!parObj.getparfloat("vmin", &vmin)) vmin = 5000;
	if (!parObj.getparfloat("padd", &padd)) padd = 25.0;
	padd = 1.0+padd/100.0;
	
	/* Get "key" values */
	nd=parObj.countparval("key");
	parObj.getparstringarray("key",key);

	/* get types and indexes corresponding to the keys */
	{ int ikey;
		for (ikey=0; ikey<nd; ++ikey) {
			type[ikey]=hdtype(key[ikey]);
			index[ikey]=getindex(key[ikey]);
		}
	}

	dx = ealloc1float(nd);
	CSMUSTGETPARFLOAT("dx",(float *)dx);
	
	if (!parObj.getparfloat("fac", &fac)) fac = 1.0;
	fac = MAX(fac,1.0);

	/* get the first record */
	rec_o = get_gather(&keyg,&typeg,&valg,&nt,&ntr,&dt,&first);
	if(ntr==0) throw cseis_geolib::csException("Can't get first record\n");
	
	/* set up the fft */
	nfft = npfar(nt*padd);
	if (nfft >= SU_NFLTS || nfft >= PFA_MAX)
		 	throw cseis_geolib::csException("Padded nt=%d--too big", nfft);
	nf = nfft/2 + 1;
	snfft=1.0/nfft;
	d1 = 1.0/(nfft*dt);
	
	rt = ealloc1float(nfft);
	ctmix = ealloc1complex(nf);
	
	
	do {
		ng++;
		 	
		fd = ealloc2complex(nf,ntr); 
		memset( (void *) ctmix, (int) '\0', nf*sizeof(complex));
		
		itrm = ealloc1float(ntr);
		rtrm = ealloc1float(ntr);
		wht = ealloc1float(ntr);

		/* transform the data into FX domain */
		{ unsigned int itr;
			for(itr=0;itr<ntr;itr++) {
				memcpy( (void *) rt, (const void *) (*rec_o[itr]).data,nt*FSIZE);
				memset( (void *) &rt[nt], (int) '\0', (nfft - nt)*FSIZE);
				pfarc(1, nfft, rt, fd[itr]);
			
			}
		}
		
		/* Do the mixing */
		{ unsigned int imx=0,itr,ifr;
		  float dist;
		  
		  	
			/* Find the trace to mix */
			for(itr=0;itr<ntr;itr++) 
				if((*rec_o[itr]).mark) {
					imx = itr;
					break;
				}
			
			memcpy( (void *) ctmix, (const void *) fd[imx],nf*sizeof(complex));
			
			/* Save the header */
			memcpy( (void *) &tr, (const void *) rec_o[imx],HDRBYTES);
 		  	
			/* weights */
			wht[imx] = 1.0;
			for(itr=0;itr<imx;itr++) {
				 dist=n_distance(rec_o,index,type,dx,nd,imx,itr);
				 wht[itr] = MIN(1.0/dist,1.0);
				 wht[itr] = 1.0;
			}
			
			for(itr=imx+1;itr<ntr;itr++) {
				 dist=n_distance(rec_o,index,type,dx,nd,imx,itr);
				 wht[itr] = MIN(1.0/dist,1.0);
				 wht[itr] = 1.0;
			}
				 
			
			/* Do the alpha trim for each trace */			
			for(ifr=0;ifr<nf;ifr++) {
 		  		for(itr=0;itr<ntr;itr++) {
					itrm[itr] = fd[itr][ifr].i;
					rtrm[itr] = fd[itr][ifr].r;
				}
				ctmix[ifr].i = alpha_trim_w(itrm,wht,ntr,trimp);
				ctmix[ifr].r = alpha_trim_w(rtrm,wht,ntr,trimp);
			}
			
					
		}
		
		
		{ unsigned int it;
			pfacr(-1, nfft, ctmix, rt);
				for(it=0;it<nt;it++) 		
					tr.data[it]=rt[it]*snfft;
		}
			
		free2complex(fd);

		{ unsigned int itr;
			for(itr=0;itr<ntr;itr++) {
				free1((void *)rec_o[itr]);
			}
		}
		
		su2cs->putTrace(&tr);
		
	    	rec_o = get_gather(&keyg,&typeg,&valg,&nt,&ntr,&dt,&first);
		
		fprintf(stderr," %d %d\n",ng,ntr);
		
		free1float(rtrm);
		free1float(itrm);
		free1float(wht);
		
	} while(ntr);
		
	
	free1float(rt);

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


float alpha_trim(float *a,int n,float p)
/* Alpha trim of array a
   the value a[irp] is replaced by the mean of the 
   a values, but the largers and smallest p% of a 
   values are left out fo mean computation 
   ,array a is destroyed */
   
{
	float mean=0;
	unsigned int ist=NINT(n*(p/100.0));
	unsigned int i,ic;
	
	qksort(n,a);
	
	for(i=ist,ic=0;i<n-ist;i++,ic++) 
		mean += a[i];
	
	return(mean/(float)ic);
}

float alpha_trim_w(float *a,float *w,int n,float p)
/******************************************************************
alpha_trim_w - Weighted Alpha trim of array 
******************************************************************
Notes:
 the value a[irp] is replaced by the mean of the 
 a values, but the largers and smallest p% of a 
 values are left out fo mean computation,array a is destroyed 
******************************************************************
Author: Balasz Nemeth, given to CWP in 2008 by Potash Corporation
******************************************************************/
{
	float mean=0;
	unsigned int ist=NINT(n*(p/100.0));
	unsigned int i,ic;
	int *ind=NULL;
	
	ind = ealloc1(n,sizeof(unsigned int));
	
	for(i=0;i<n;++i) ind[i]=i;
	
	qkisort(n,a,ind);
	
	for(i=ist,ic=0;i<n-ist;i++,ic++) 
		mean += a[ind[i]]*w[ind[i]];
		
	free1(ind);
	
	return(mean/(float)ic);
}




float n_distance(segy **rec_o,int *index,cwp_String *type,float *dx,unsigned int nd,unsigned int imx,unsigned int itr)
{

	float du;
	float o,x;
	int i;
	double dist=0;
	Value val;
	
	for(i=0;i<nd;i++) {
 		gethval(rec_o[imx], index[i],&val);
		o = vtof(type[i],val);
 		gethval(rec_o[itr], index[i],&val);
		x = vtof(type[i],val);
		du = (o-x)/dx[i];
		dist += du*du;
	}
	return((float)sqrt(dist));	
}
	

} // END namespace
