/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUMIXGATHER:  $Date: June 2000  */

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
std::string sdoc_sumixgathers =
"									"
" SUMIXGATHERS - mix two gathers					"
" 									"
" sumixgathers file1 file2 > stdout [optional parameters]		"
" 									"
" Required Parameters:							"
" ntr=tr.ntr	if ntr header field is not set, then ntr is mandatory	"
"									"
" Optional Parameters: none						"
"									"
" Notes: Both files have to be sorted by offset				"
" Mixes two gathers keeping only the traces of the first file		"
" if the offset is the same. The purpose is to substitute only		"
" traces non existing in file1 by traces interpolated store in file2. 	" 
"									"
" Example. If file1 is original data file and file 2 is obtained by	"
" resampling with Radon transform, then the output contains original  	"
" traces with gaps filled						"
" 									"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sumixgathers {


/* Credits:
 *	Daniel Trad. UBC
 * Trace header fields accessed: ns, dt, ntr
 * Copyright (c) University of British Columbia, 1999.
 * All rights reserved.
 */

/**************** end self doc ***********************************/

segy tr,tr2; 
void* main_sumixgathers( void* args )
{
	int j,ih;
	register int it;
	int nt;
	int ntr;
	int nh;
	int flag;
	float *h;
	float scale=1;
	int scaling;

	FILE *fp1=NULL;		/* file pointer for first file		*/
	FILE *fp2=NULL;		/* file pointer for second file		*/
	FILE *tracefp=NULL;	/* fp for trace storage file		*/
	FILE *headerfp=NULL;	/* fp for header storage file		*/

	/* Initialize */

	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sumixgathers );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Open two files given as arguments for reading */
	fp1 = efopen(argv[1], "r");
	fp2 = efopen(argv[2], "r");  
	tracefp = etmpfile();
	headerfp = etmpfile();
	if (!parObj.getparint("scaling",&scaling)) scaling=0;  
	if (!fgettr(fp1,&tr)) throw cseis_geolib::csException("can't read first trace");

	nt = (int) tr.ns;  
	if (!parObj.getparint("ntr",&ntr)) ntr=tr.ntr;  
        parObj.checkpars();
	if (!ntr)
		throw cseis_geolib::csException("ntr neither set in header nor parObj.getparred!");
	h=ealloc1float(ntr);

	j=0;
	do{	/* Loop over traces */	
		efwrite(&tr,HDRBYTES,1,headerfp);
		efwrite(tr.data,FSIZE, nt, tracefp);  		
		h[j]=tr.offset;
		j++;
	} while(fgettr(fp1, &tr));

	erewind(tracefp);
	erewind(headerfp);
	nh=j;

	warn("nh=%d",nh);

	scale=1;

	while(fgettr(fp2, &tr2)){
		flag=1;
		for (ih=0;ih<nh;ih++){
			if (h[ih]>=0) 
			if ((tr2.offset>0.999*h[ih])&&(tr2.offset<1.001*h[ih])) flag=0;
			if (h[ih]<0) 
			if ((tr2.offset<0.99*h[ih])&&(tr2.offset>1.01*h[ih])) flag=0;
		}
		if (flag==1){
			if (scaling && fabs(tr2.offset) > 0){
		  	  scale=(1+0.03*(fabs(tr2.offset)/1000.));
			  fprintf(stderr,"tr2.offset=%d,scale=%f\n",tr2.offset,scale);
			  for (it=0;it<nt;it++) tr2.data[it]*=scale;
			}
		su2cs->putTrace(&tr2);
		}
	}
	
	for (ih=0;ih<nh;ih++){
		efread(&tr,HDRBYTES,1,headerfp);
		efread(tr.data,FSIZE, nt, tracefp);
		su2cs->putTrace(&tr);
	}
	

	efclose(fp1);
	efclose(fp2);
	efclose(tracefp);
	efclose(headerfp);
	free1float(h);

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
