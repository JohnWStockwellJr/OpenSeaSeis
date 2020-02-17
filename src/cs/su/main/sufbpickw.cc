/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUFBPICKW: $Revision: 1.3 $ ; $Date: 2011/11/16 23:13:27 $  */

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

/*********************** self documentation *****************************/
std::string sdoc_sufbpickw =
" SUFBPICKW - First break auto picker				"
"								"
"   sufbpickw < infile >outfile					"
"								"
" Required parameters:						"
"  none								"
" Optional parameters:						"
" keyg=ep						 	"
" window=.03	Length of forward and backward windows (s)	"
" test=1	Output the characteristic function	 	"
"		This can be used for testing window size	"
" Template							"
" o=		offset...				  	"
" t=		time pairs for defining first break search	"
"			window centre				"
" tdv=.05	Half length of the search window		"
"								"
" If the template is specified the maximum value of the		"
" characteristic function is searched in the window		"
"  defined by the template only.Default is the whole trace.	"
"							  	"
"  The time of the pick is stored in header word unscale 	"
"								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sufbpickw {

   
/* segy data  */
segy *trp;				/* SEGY trace array */
segy trtp;				/* SEGY trace */
void* main_sufbpickw( void* args )
{
	
	segy **rec_o;	   /* trace header+data matrix */  
	
	cwp_String keyg;
	cwp_String typeg;		
	Value valg;
		   	
	int first=0;		/* true when we passed the first gather */
	int ng=0;
	float dt;
	int nt;
	int ntr;

	unsigned int np;	/* Number of points in pick template */
	float *t=NULL;		/* array defining pick template times */
	float *o=NULL;		/* array defining pick template offsets */
	
	float window;
	int iwindow;
	int *itimes;
	int *itimes2;
	float *offset;
	float tdv;
	int itdv;
	float *find;
	int nowindow=0;
	int test=0;
	
	FILE *ttp;
		
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sufbpickw );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	
	if (!parObj.countparval("t")) {
		np=2;
		nowindow=1;
	} else {
		np=parObj.countparval("t");
	}
	
	if(!nowindow) { 
		t  = ealloc1float(np);
		o  = ealloc1float(np);
	
		if( np == parObj.countparval("o")) {
			parObj.getparfloat("t",t);
			parObj.getparfloat("o",o);
		} else {
			throw cseis_geolib::csException(" t and o has different number of elements\n");
		}
	}



	ttp = efopen("test.su","w");
		
	if (!parObj.getparstring("keyg", &keyg)) keyg ="ep";
	
	if (!parObj.getparfloat("window",&window)) window =0.02;
	if (!parObj.getparfloat("tdv",&tdv)) tdv = -1.0;
	if (!parObj.getparint("test",&test)) test = 1;
	
        parObj.checkpars();
	/* get information from the first header */
	rec_o = get_gather(&keyg,&typeg,&valg,&nt,&ntr,&dt,&first);
	
	iwindow=NINT(window/dt);
	if(tdv==-1.0) {
		itdv=nt;
	} else {
		itdv=NINT(2.0*tdv/dt);
	}
	
	if(ntr==0) throw cseis_geolib::csException("Can't get first record\n");
	do {
		ng++;
		
		itimes = ealloc1int(ntr);
		itimes2 = ealloc1int(ntr);
		offset = ealloc1float(ntr);
		
		/* Phase 1 */
		/* Loop through traces */
		{ int itr,ifbt;
		  int it;
		  float fbt;
		  float ampb;		
		  float ampf;
		  float *wf,*wb;
		  float *cf;		/* Characteristic function */
		  
		  
		  cf=ealloc1float(nt);
		  find=ealloc1float(nt);
		  for(it=0;it<nt;it++)
		  	find[it]=(float)it;
		  
		  for(itr=0;itr<ntr;itr++) {
		  
		  	memset( (void *) cf, (int) '\0', nt*FSIZE);
			
			if(nowindow) {
				ifbt=0;
			} else {
		  	/* Linear inperpolation of estimtated fb time */
				offset[itr] =(*rec_o[itr]).offset; 
		  		intlin(np,o,t,t[0],t[np-1],1,&offset[itr],&fbt);
		  		ifbt = NINT(fbt/dt);
			}
			
		  
			wb = &(*rec_o[itr]).data[0];
			wf = &(*rec_o[itr]).data[iwindow];
			ampb = sasum(iwindow,wb,1)+FLT_EPSILON;
			ampf = sasum(iwindow,wf,1)+FLT_EPSILON;
		  	for(it=iwindow;it<nt-iwindow-1;it++) {
				cf[it] = ampf/ampb;
				/* setup next window */
				ampb -= fabs((*rec_o[itr]).data[it-iwindow]);
				ampf -= fabs((*rec_o[itr]).data[it]);
				ampb += fabs((*rec_o[itr]).data[it]);
				ampf += fabs((*rec_o[itr]).data[it+iwindow]);
			}
			/* Smooth the characteristic function */
			smooth_segmented_array(&find[iwindow],&cf[iwindow],nt-iwindow-1,iwindow,1,5);
		  	
			/* find the maximum*/
			it=MIN(MAX(ifbt-itdv/2,0),nt-1);
			

			wb = &cf[it];
			
			itimes[itr] = isamax(itdv,wb,1)+it;
			
			/* Final check */
			if(itimes[itr] == ifbt-itdv/2 || itimes[itr] == ifbt+itdv/2)
				itimes[itr]=0;
				 
			(*rec_o[itr]).unscale=dt*itimes[itr];
			
			if(test)
				memcpy( (void *) &(*rec_o[itr]).data[0], (const void *) cf, nt*FSIZE);

		  }
		  free1float(cf);
		  free1float(find);
		}

		free1int(itimes);
		free1int(itimes2);
		free1float(offset);
		
		rec_o = put_gather(rec_o,&nt,&ntr);
		rec_o = get_gather(&keyg,&typeg,&valg,&nt,&ntr,&dt,&first);
	} while(ntr);
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
