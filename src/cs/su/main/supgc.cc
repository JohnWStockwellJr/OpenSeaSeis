/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

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
std::string sdoc_supgc =
" SUPGC   -   Programmed Gain Control--apply agc like function	"
"              but the same function to all traces preserving	"
"              relative amplitudes spatially.			"
" Required parameter:						"
" file=             name of input file				"
"								"
" Optional parameters:						"
" ntrscan=200       number of traces to scan for gain function	"
" lwindow=1.0       length of time window in seconds		"
"								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace supgc {


/*
 * Author: John Anderson (visitor to CWP from Mobil)
 *
 * Trace header fields accessed: ns, dt
 */

/**************** end self doc ********************************/

segy tr;	/* Input and output trace data of length nt */
void* main_supgc( void* args )
{
	FILE *fp;
	char *file;
	int ntrscan;
	int icount,j,k,kk,nt,lw;
	float lwindow,dt,sum;
	float *gain,*g;

	/* hook up getpars */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_supgc );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* get parameters */
	if (!parObj.getparint("ntrscan",&ntrscan)) ntrscan = 200;
        if(!parObj.getparstring("file",&file)) file=NULL;
	if(!parObj.getparfloat("lwindow",&lwindow)) lwindow=1.0;

        parObj.checkpars();

	/* get info from first trace */
        if( (fp = fopen(file,"r"))==NULL) 
		throw cseis_geolib::csException("could not open input file");
        if(!fgettr(fp,&tr)) throw cseis_geolib::csException("Can't get first trace");
	nt=tr.ns;
	dt = ((double) tr.dt)/1000000.0;
	lw=NINT(0.5*lwindow/dt);
	warn("window length= %d samples",lw); 

	gain=ealloc1float(nt);
	g=ealloc1float(nt);

	for(j=0;j<nt;j++) gain[j]=0.;	
	icount=0;
	do{
		for(j=0;j<nt;j++)  gain[j]+=fabs(tr.data[j]);
		icount++;
		
	} while(fgettr(fp,&tr) && icount<ntrscan);

	rewind(fp);

	for(j=0;j<nt;j++) {
		sum=0.;
		kk=0;
		for(k=MAX(j-lw,0);k<MIN(j+lw,nt);k++){
			kk++;
			sum+=gain[k];
		}
		if(sum==0.) sum=1.;
		g[j]=kk*icount/sum;
	}
	
	warn("scan covered %d traces",icount);		

	if(!fgettr(fp,&tr))throw cseis_geolib::csException("failure after rewind");
	do {
		for(j=0;j<nt;j++) tr.data[j]*=g[j];

		su2cs->putTrace(&tr);

	} while(fgettr(fp,&tr));


	fclose(fp);
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
