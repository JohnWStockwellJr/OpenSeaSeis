/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUADDSTATICS: $Revision: 1.6 $ ; $Date: 2011/11/16 23:16:23 $	*/

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
#include "cwp.h"
#include "header.h"
#include "segy.h"

/*********************** self documentation **********************/
std::string sdoc_suaddstatics =
"									"
" SUADDSTATICS - ADD random STATICS on seismic data			"
"								 	"
" suaddstatics required parameters [optional parameters] > stdout	"
"									"
" Required parameters:							"
" shift=		the static shift will be generated 	 	"
"			randomly in the interval [+shift,-shif] (ms)	"
" sources=		number of source locations			"
" receivers=		number of receiver locations			"
" cmps=			number of common mid point locations		"
" maxfold=		maximum fold of input data			"
" datafile=		name and COMPLETE path of the input file	"
"									"
" Optional parameters:							"
" dt=tr.dt			time sampling interval (ms)		"
" seed=getpid()		 seed for random number generator		"
" verbose=0			=1 print useful information		"
"									"
" Notes:								"
" Input data should be sorted into cdp gathers.				"
"									"
" SUADDSTATICS applies static time shifts in a surface consistent way on"
" seismic data sets. SUADDSTATICS writes the static time shifts in the  "
" header field TSTAT. To perform the actual shifts the user should use 	"
" the program SUSTATIC after SUADDSTATICS. SUADDSTATICS outputs the	"
" corrupted data set to stdout.						"
"									"
" Header field used by SUADDSTATICS: cdp, sx, gx, tstat, dt.		"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suaddstatics {


/*
 * Credits: CWP Wences Gouveia, 11/07/94,  Colorado School of Mines
 */
/************************ end self doc ***********************************/

segy tr;		/* SEGY DATA */

extern void sranuni();
extern double Randdouble();

void* main_suaddstatics( void* args )
{
	double max_stat,min_stat; /* limits of statics shift	*/

	float *vect=NULL;	/* floating point representation*/
	float ***data=NULL;	/* seismic data			*/
	float dt=0;		/* time sampling		*/
	float trace_static;	/* total statics, in SECONDS	*/
	float receiver_statics;
	float source_statics;   /* statics components		*/
	float shift;		/* shift to be applied 		*/

	int ***pointer=NULL;	/* set of pointers		*/
	int *source=NULL;
	int *receiver=NULL;	/* """ "" """"""""		*/
	int *foldcmp=NULL;	/* fold of each CMP		*/

	int isource=0;		/* source counter		*/
	int itrace=0; 		/* trace counter		*/
	int ireceiver=0;	/* receiver counter 		*/
	int icmp=0;		/* cmp counter			*/
	int cmp_current=0;	/* current cmp			*/
	int cmp_last=0;		/* last cmp			*/

	int seed;		/* seed for random number generator */
	int i;			/* counters and similar stuff	*/
	int ncmp;		/* number of cmps		*/
	int nsources;		/* number of sources		*/
	int nreceivers;		/* number of receivers		*/
	int maxfold;		/* maximum fold of the data	*/
	int nsamples;		/* number of samples / trace	*/
	int ntraces;		/* used in reading input data	*/
	int verbose;		/* dialogue			*/

	char *datafile=NULL;		/* input file			*/
	FILE *fp=NULL, *hfp=NULL;	/* files ....			 */

	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suaddstatics );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get parameters */
	/* Required parameters */
	if (!parObj.getparint("sources", &nsources)) 
		throw cseis_geolib::csException("Specify number of sources!");
	if (!parObj.getparint("receivers", &nreceivers))
		throw cseis_geolib::csException("Specify number of receivers!");
	if (!parObj.getparint("cmps", &ncmp))
		throw cseis_geolib::csException("Specify number of cmps!");
	if (!parObj.getparint("maxfold", &maxfold))
		throw cseis_geolib::csException("Specify maximum fold of the data!");
	if (!parObj.getparstring("datafile", &datafile))
		throw cseis_geolib::csException("Specify data file name!");

	/* Optional Parameters */
	if (!parObj.getparint("seed", &seed)) seed = getpid();
	if (!parObj.getparint("verbose", &verbose)) verbose = 0;
	if (!parObj.getparfloat("dt", &dt))  dt = 0.;
	if (!parObj.getparfloat("shift", &shift))
		throw cseis_geolib::csException("Specify time shift (sec)!");
	if (verbose) warn("Reading input data set");
        parObj.checkpars();

	/* Open data file */
	fp = fopen(datafile,"r");
	if ((fp = fopen(datafile, "r")) == NULL)
		throw cseis_geolib::csException("SUADDSTATICS: can't open input data file");

	if (!fgettr(fp,&tr)) throw cseis_geolib::csException("can't get first trace");
	nsamples = tr.ns;
	if (dt == 0.) dt = (float) tr.dt / 1000.0;
	rewind(fp);

	/* Allocate memory */
	pointer = alloc3int(2,maxfold,ncmp);
	vect = alloc1float(nsources + nreceivers);
	source = alloc1int(nsources);
	receiver = alloc1int(nreceivers);
	foldcmp = alloc1int(ncmp);
	data = alloc3float(nsamples, maxfold, ncmp);

/*
	Reading data and headers. 
	And copying traces and headers to temporary files 
	NOTE THAT CMPs MUST BE IN ASCENDING ORDER 
*/
	hfp = tmpfile();
	ntraces = 0;
	isource = 0;
	ireceiver = 0;
	icmp = -1;
	cmp_current = 0;
	cmp_last = 0;
	while (fgettr(fp,&tr)) { 
	 	cmp_current = tr.cdp;
		if (cmp_current != cmp_last)
		{
			if (icmp >= 0) foldcmp[icmp] = itrace;
			icmp++;
			itrace = 0;
			cmp_last = cmp_current;
		}	
	 	pointer[icmp][itrace][0] = tr.sx;	
		pointer[icmp][itrace][1] = tr.gx;

		if (icmp == 0 && itrace == 0)
			source[isource] = pointer[icmp][itrace][0];

		else if (!(icmp == 0 && itrace == 0) && source[isource] != pointer[icmp][itrace][0])
		{
			i = isource - 1;
			while (i >= 0 && source[i] != pointer[icmp][itrace][0])
			i--;
			if (i < 0)
			{
				isource++;
                                if(isource >= nsources){
                                   break;
                                }
				source[isource] = pointer[icmp][itrace][0];
			}
		}

		if (icmp == 0 && itrace == 0)
			receiver[ireceiver] = pointer[icmp][itrace][1];

		else if (!(icmp == 0 && itrace == 0) && receiver[ireceiver] != pointer[icmp][itrace][1]) {
			i = ireceiver - 1;
			while (i >= 0 && receiver[i] != pointer[icmp][itrace][1])
			i--;
			if (i < 0) {
				ireceiver++;
                                if( ireceiver >= nreceivers){
                                   break;
                                }
				receiver[ireceiver] = pointer[icmp][itrace][1];
			}
		}
		efwrite(&tr,HDRBYTES,1,hfp);
		memcpy(data[icmp][itrace],tr.data,nsamples*FSIZE);
		itrace++;
		ntraces++;
	};
	foldcmp[icmp] = itrace;	/* last CMP */
	/* temporal sampling */
	erewind(hfp); 

	if (verbose)
	   fprintf(stderr,"%d TRACES were read in %d CMPs\n",ntraces,icmp+1);
/*
	Enhancing the pointers
*/
	for (icmp = 0; icmp < ncmp; icmp++) {
		for (itrace = 0; itrace < foldcmp[icmp]; itrace++) {
			isource = 0;
			while (source[isource] != pointer[icmp][itrace][0])
				isource++;
			pointer[icmp][itrace][0] = isource;

			ireceiver = 0;
			while (receiver[ireceiver] != pointer[icmp][itrace][1])
				ireceiver++;	
			pointer[icmp][itrace][1] = ireceiver;
		}
	}

	free1int(source);
	free1int(receiver);	/* Neither for this pointers	   */

	sranuni(seed);
	min_stat = -shift / dt;
	max_stat = shift / dt;

	for (i = 0; i < nsources + nreceivers; i++)
		vect[i] = Randdouble(min_stat,max_stat);

/*
	And reading and writing the data and headers with the field 
	tstat properly filled
*/
	for (icmp = 0; icmp < ncmp; icmp++) {
		if (verbose)
			fprintf(stderr,"Applying fake statics to CMP# %d\n",icmp);

		for (itrace = 0; itrace < foldcmp[icmp]; itrace++) {
			efread(&tr,HDRBYTES,1,hfp);
/*
	trace_static must be in MILISECONDS
*/
			source_statics = vect[pointer[icmp][itrace][0]] * dt;
			receiver_statics = vect[pointer[icmp][itrace][1] + nsources] * dt;
			trace_static = source_statics + receiver_statics; 
			tr.tstat = NINT(trace_static); 
			memcpy(tr.data,data[icmp][itrace],nsamples*FSIZE);
			su2cs->putTrace(&tr);
		}
	}
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

double Randdouble(dlow,dhigh)
double dlow, dhigh;
{
	float random;
	double drandom;
	
	random = franuni();
	drandom = (double) (random * (dhigh - dlow) + dlow);
	return (drandom);
}
/*** end of file ***/

} // END namespace
