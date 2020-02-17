/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUFWMIX: $Revision: 1.4 $ ; $Date: 2011/11/16 23:09:52 $		*/

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
std::string sdoc_sufwmix =
" SUFWMIX -  FX domain multidimensional Weighted Mix			"
"									"
"	sufwmix < stdin > stdout [optional parameters]			"
"									"
" Required parameters:							"
" key=key1,key2,..	Header words defining mixing dimension		"
" dx=d1,d2,..		Distance units for each header word		"
" Optional parameters:							"
" keyg=ep		Header word indicating the start of gather	"
" vf=0			=1 Do a frequency dependent mix			"
" vmin=5000		Velocity of the reflection slope		"
"			than should not be attenuated			"
" Notes:								"
" Trace with the header word mark set to one will be			"
" the output trace 							"
"  (a work in progress)							"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sufwmix {


/*
 * Credits:  
 *
 *  Potash Corporation: Balazs Nemeth, Saskatoon Saskatchewan CA,
 *   given to CWP in 2008
 *
 */

/**************** end self doc ********************************/


#define PFA_MAX 720720  /* Largest allowed nfft	   */
#define PIP2 1.570796327
#define AMPSP(c) rcabs(c)
#define PHSSP(c) atan2(c.i,c.r)

segy tr;

/* function prototype of subroutine used internally */
float n_distance(segy **rec_o,int *index, cwp_String *type,
			float *dx,unsigned int nd,
			unsigned int imx,unsigned int itr);

int verbose;

void* main_sufwmix( void* args )
{
	cwp_String keyg;	/* header key word from segy.h		*/
	cwp_String typeg;	/* ... its type				*/
	Value valg;
	cwp_String key[SU_NKEYS];	/* array of keywords		 */
	cwp_String type[SU_NKEYS];	/* array of keywords		 */
	int index[SU_NKEYS];		/* name of type of getparred key  */
	
	
	int first=0;		/* true when we passed the first gather */
	int ng=0;
	float dt;
	int nt;		
	int ntr;		/* number of traces in a gather */
	
	int nfft=0;		/* length of padded array */
	float snfft;		/* scale factor for inverse fft */
	int nf=0;		/* number of frequencies */
	float d1;		/* frequency sampling int. */
	float *rt;		/* real trace */
	complex *ctmix;		/* complex trace */
	complex **fd;		/* frequency domain data */

	
	float padd;
	
	int nd;			/* number of dimensions */
	float *dx;
	float fac;
	float vmin;
	int vf;
		
	segy **rec_o;		/* trace header+data matrix */	

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sufwmix );
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

        parObj.checkpars();

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
		  float weight=0,dist=0,weightsum=0,dlim;
		  
		  	
			/* Find the trace to mix */
			for(itr=0;itr<ntr;itr++) 
				if((*rec_o[itr]).mark) {
					imx = itr;
					break;
				}
			
			memcpy( (void *) ctmix, (const void *) fd[imx],nf*sizeof(complex));
			
			/* Save the header */
			memcpy( (void *) &tr, (const void *) rec_o[imx],HDRBYTES);
				
			weightsum=0;
				
			/* Do the mixing up to the imx trace */			
		  	for(itr=0;itr<imx;itr++) {
				
				/* mixing weight 1.0/distance^2 */
				dist=n_distance(rec_o,index,type,dx,nd,imx,itr);
				weight  = 1.0/(dist/fac+1.0);
				weightsum=+weight;
				
				for(ifr=0;ifr<nf;ifr++) {
					if(vf) {
						dlim=vmin/(2.0*ifr*d1);
						if(dlim>dist) {
							ctmix[ifr].r += weight*fd[itr][ifr].r;
							ctmix[ifr].i += weight*fd[itr][ifr].i;
						}
					} else {
						ctmix[ifr].r += weight*fd[itr][ifr].r;
						ctmix[ifr].i += weight*fd[itr][ifr].i;
					}
				}
			}
				
			/* Do the mixing from imx till the end  */			
		  	for(itr=imx+1;itr<ntr;itr++) {
				
				/* mixing weight 1.0/distance^2 */
				dist=n_distance(rec_o,index,type,dx,nd,imx,itr);
				weight  = 1.0/(dist/fac+1.0);
				weightsum=+weight;
				
				for(ifr=0;ifr<nf;ifr++) {
					if(vf) {
						dlim=vmin/(2.0*ifr*d1);
						if(dlim>dist) {
							ctmix[ifr].r += weight*fd[itr][ifr].r;
							ctmix[ifr].i += weight*fd[itr][ifr].i;
						}
					} else {
						ctmix[ifr].r += weight*fd[itr][ifr].r;
						ctmix[ifr].i += weight*fd[itr][ifr].i;
					}
				}
			}
			
					
			/* Normalize */
			for(ifr=0;ifr<nf;ifr++) 
				crmul(ctmix[ifr],1.0/weightsum);
					
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
