/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* FX domain filter corelation coefficient filter*/


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


#define PFA_MAX 720720  /* Largest allowed nfft           */
#define PIP2 1.570796327
#define AMPSP(c) rcabs(c)
#define PHSSP(c) atan2(c.i,c.r)


/*********************** self documentation **********************/
std::string sdoc_succfilt =
" SUCCFILT -  FX domain Correlation Coefficient FILTER			"
"                       						"
"   sucff < stdin > stdout [optional parameters]			"
"                       						"
" Optional parameters:							"
" cch=1.0		Correlation coefficient high pass value		"
" ccl=0.3		Correlation coefficient low pass value		"
" key=ep		ensemble identifier				"
" padd=25		FFT padding in percentage			"
"                       						"
" Notes:                       						"
" This program uses \"get_gather\" and \"put_gather\" so requires that	"
" the  data be sorted into ensembles designated by \"key\" with the ntr"
" field set to the number of traces in each respective ensemble.  	"
"                       						"
" Example:                     						"
" susort ep offset < data.su > datasorted.su				"
" suputgthr dir=Data verbose=1 < datasorted.su				"
" sugetgthr dir=Data verbose=1 > dataupdated.su				"
" succfilt  < dataupdated.su > ccfiltdata.su				"
"                       						"
" (Work in progress, editing required)                 			"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace succfilt {


/* 
 * Credits:
 *  Potash Corporation: Balazs Nemeth, Saskatoon Canada. c. 2008
 */


/**************** end self doc ********************************/

segy tr;
segy tr2;

void* main_succfilt( void* args )
{
	cwp_String key;		/* header key word from segy.h		*/
	cwp_String type;	/* ... its type				*/
	Value val;

	segy **rec_o=NULL;	/* trace header+data matrix */

	int first=0;		/* true when we passed the first gather */
	int ng=0;		/* counter of gathers */
	float dt;
	int nt;
	int ntr;
	
	int nfft=0;		/* lenghth of padded array */
	float snfft;		/* scale factor for inverse fft */
	int nf=0;		/* number of frequencies */
	float *rt=NULL;		/* real trace */
	complex **fd=NULL;	/* frequency domain data */
	float **cc=NULL;	/* correlation coefficinet matrix */
	
	float padd;
	float cch;
	float ccl;

	int verbose=0;
		
	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_succfilt );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	
	if (!parObj.getparstring("key", &key))	key = "ep";
	if (!parObj.getparfloat("padd", &padd)) padd = 25.0;
	padd = 1.0+padd/100.0;
	
	if (!parObj.getparfloat("cch", &cch)) cch = 1.0;
	if (!parObj.getparfloat("ccl", &ccl)) ccl = 0.3;
	if (!parObj.getparint("verbose",&verbose)) verbose = 0;
	
	/* get the first record */
	rec_o = get_gather(&key,&type,&val,&nt,&ntr,&dt,&first);
	if(ntr==0) throw cseis_geolib::csException("Can't get first record\n");
	
	/* set up the fft */
	nfft = npfar(nt*padd);
        if (nfft >= SU_NFLTS || nfft >= PFA_MAX)
               	throw cseis_geolib::csException("Padded nt=%d--too big", nfft);
        nf = nfft/2 + 1;
        snfft=1.0/nfft;
	
	rt = ealloc1float(nfft);
	
	do {
		ng++;

		fd = ealloc2complex(nf,ntr); 
       		cc = ealloc2float(nf,ntr);

		/* transform the data into FX domain */
		{ unsigned int itr;
			for(itr=0;itr<ntr;itr++) {
				memcpy( (void *) rt, (const void *) (*rec_o[itr]).data,nt*FSIZE);
                		memset( (void *) &rt[nt],0, (nfft - nt)*FSIZE);
				
				pfarc(1, nfft, rt, fd[itr]);
			}
		}
		
		/* Compute correlation coefficients */
		{ unsigned int itr,ifr;
			for(itr=0;itr<ntr-1;itr++) {
				for(ifr=0;ifr<nf-1;ifr++) { 
					cc[itr][ifr] = cos(PHSSP(fd[itr][ifr])-PHSSP(fd[itr+1][ifr])); 
				}			
			}
		
		}
		
		/* Filter */
		{ unsigned int itr,ifr;
			for(itr=0;itr<ntr-1;itr++) {
				for(ifr=0;ifr<nf-1;ifr++) { 
					if(cc[itr][ifr]> cch || cc[itr][ifr]<ccl) {
						fd[itr][ifr].r = 0.0; 
						fd[itr][ifr].i = 0.0;
					} 
				}			
			}
		
		}
		
		{ unsigned int itr,it;
			for(itr=0;itr<ntr;itr++) {
				
				pfacr(-1, nfft, fd[itr], rt);
				
				for(it=0;it<nt;it++) 		
                			(*rec_o[itr]).data[it]=rt[it]*snfft;
			}
		}
			
		free2complex(fd);
		free2float(cc);

	    	rec_o = put_gather(rec_o,&nt,&ntr);
	    	rec_o = get_gather(&key,&type,&val,&nt,&ntr,&dt,&first);
		
		if (verbose) warn(" %d %d\n",ng,ntr);
		
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


} // END namespace
