/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUADDNOISE: $Revision: 1.48 $ ; $Date: 2011/11/12 00:22:43 $		*/

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
extern "C" {
  #include <time.h>

}
extern "C" {
  #include <signal.h>

}

/*********************** self documentation ******************************/
std::string sdoc_suaddnoise =
" 									"
" SUADDNOISE - add noise to traces					"
" 									"
" suaddnoise <stdin >stdout  sn=20  noise=gauss  seed=from_clock	"
" 									"
" Required parameters:							"
" 	if any of f=f1,f2,... and amp=a1,a2,... are specified by the user"
"	and if dt is not set in header, then dt is mandatory		"
" 									"
" Optional parameters:							"
" 	sn=20			signal to noise ratio			"
" 	noise=gauss		noise probability distribution		"
" 				=flat for uniform; default Gaussian	"
" 	seed=from_clock		random number seed (integer)		"
"	f=f1,f2,...		array of filter frequencies (as in sufilter)"
"	amps=a1,a2,...		array of filter amplitudes		"
" 	dt= (from header)	time sampling interval (sec)		"
"	verbose=0		=1 for echoing useful information	"
" 									"
" 	tmpdir=	 if non-empty, use the value as a directory path	"
"		 prefix for storing temporary files; else if the	"
"	         the CWP_TMPDIR environment variable is set use		"
"	         its value for the path; else use tmpfile()		"
" 									"
" Notes:								"
" Output = Signal +  scale * Noise					"
" 									"
" scale = (1/sn) * (absmax_signal/sqrt(2))/sqrt(energy_per_sample)	"
" 									"
" If the signal is already band-limited, f=f1,f2,... and amps=a1,a2,...	"
" can be used, as in sufilter, to bandlimit the noise traces to match	"
" the signal band prior to computing the scale defined above.		"
" 									"
" Examples of noise bandlimiting:					"
" low freqency:    suaddnoise < data f=40,50 amps=1,0 | ...		"
" high freqency:   suaddnoise < data f=40,50 amps=0,1 | ...		"
" near monochromatic: suaddnoise < data f=30,40,50 amps=0,1,0 | ...	"
" with a notch:    suaddnoise < data f=30,40,50 amps=1,0,1 | ...	"
" bandlimited:     suaddnoise < data f=20,30,40,50 amps=0,1,1,0 | ...	"
" 									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suaddnoise {


/* Credits:
 *	CWP: Jack Cohen, Brian Sumner, Ken Larner
 *		John Stockwell (fixed filtered noise option)
 *
 * Notes:
 *	At S/N = 2, the strongest reflector is well delineated, so to
 *	see something 1/nth as strong as this dominant reflector
 *	requires S/N = 2*n.
 *
 * Trace header field accessed: ns
 */

/**************** end self doc *******************************************/

/* Default signal to noise ratio */
#define SN	20

/* Noise probability distributions */
#define	GAUSS	0
#define	FLAT	1


/* Prototype */
static void closefiles(void);

/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
static char bandoutfile[L_tmpnam];  /* output file for sufilter	*/
static FILE *bandoutfp;		    /* fp for output file	*/


segy tr;

void* main_suaddnoise( void* args )
{
	int nt;			/* number of points on trace		*/
	unsigned long databytes;/* ... in bytes 			*/
	int ntr;		/* number of traces			*/
	int verbose;		/* flag for echoing info		*/
	cwp_String stype;	/* noise type (gauss, flat) as string	*/
	int itype=GAUSS;	/* ... as integer (for use in switch)	*/
	float sn;		/* signal to noise ratio		*/
	time_t seed_time;       /* random number seed                   */
	unsigned int seed;	/* random number seed			*/
	int nfloats;		/* number of floats in "signal"		*/
	float *noise;		/* noise vector				*/
	float noiscale;		/* scale for noise			*/
	float absmaxsig;	/* absolute maximum in signal		*/
	float noipow;		/* a measure of noise power		*/
	cwp_String f="";	/* frequency input for sufilter		*/
	cwp_String amps="";	/* amplitude input for sufilter		*/
	char *tmpdir;		/* directory path for tmp files		*/
	cwp_Bool istmpdir=cwp_false;/* true for user given path		*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suaddnoise );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get parameters */
	if (!parObj.getparint("verbose", &verbose))	verbose = 0;

	/* Look for user-supplied tmpdir */
	if (!parObj.getparstring("tmpdir",&tmpdir) &&
	    !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
	if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
		throw cseis_geolib::csException("you can't write in %s (or it doesn't exist)", tmpdir);

		
	/* Get noise type */
	if (!parObj.getparstring("noise", &stype))	stype = "gauss";

	/* Recall itype initialized as GAUSS */
	if (STREQ(stype, "flat"))  itype = FLAT;
	else if (!STREQ(stype, "gauss"))
		throw cseis_geolib::csException("noise=\"%s\", must be gauss or flat", stype);


	/* Get signal to noise ratio */
	if (!parObj.getparfloat("sn", &sn))	sn = SN;
	if (sn <= 0) throw cseis_geolib::csException("sn=%d must be positive", sn);


	/* Set seed */
	if (!parObj.getparuint("seed", &seed)) { /* if not supplied, use clock */
	  if (((time_t)-1) == (seed_time = time((time_t *) NULL))) {
	     throw cseis_geolib::csException("time() failed to set seed");
	     seed = 0;
	  }
	  else {
	     seed = (unsigned int)seed_time;
	  }
	}
	(itype == GAUSS) ? srannor(seed) : sranuni(seed);

	/* Prepare temporary files to hold headers and data */
	if (STREQ(tmpdir,"")) {
		tracefp = etmpfile();
		headerfp = etmpfile();
		if (verbose) warn("using tmpfile() call");
	} else { /* user-supplied tmpdir */
		char directory[BUFSIZ];
		strcpy(directory, tmpdir);
		strcpy(tracefile, temporary_filename(directory));
		strcpy(headerfile, temporary_filename(directory));
		/* Trap signals so can remove temp files */
		signal(SIGINT,  (void (*) (int)) closefiles);
		signal(SIGQUIT, (void (*) (int)) closefiles);
		signal(SIGHUP,  (void (*) (int)) closefiles);
		signal(SIGTERM, (void (*) (int)) closefiles);
		tracefp = efopen(tracefile, "w+");
		headerfp = efopen(headerfile, "w+");
      		istmpdir=cwp_true;		
		if (verbose) warn("putting temporary files in %s", directory);
	}

	/* Get info from first trace */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;
	databytes = nt * FSIZE;


	/* Loop over input traces & write headers and data to tmp files */
	ntr = 0;
	do {
		++ntr;
		efwrite(&tr, 1, HDRBYTES, headerfp);
		efwrite(tr.data, 1, databytes, tracefp);
	} while (cs2su->getTrace(&tr));
	nfloats = ntr * nt;


	/* Compute absmax of signal over entire data set */
	rewind(tracefp);
	absmaxsig = 0.0;
	{ register int i;
	  for (i = 0; i < nfloats; ++i) {
		float sigval;
		efread(&sigval, FSIZE, 1, tracefp);
		absmaxsig = MAX(absmaxsig, ABS(sigval));
	  }
	}


	/* Compute noise vector elements in [-1, 1] */
	noise = ealloc1float(nfloats);
	switch (itype) {
		register int i;
	case GAUSS: /* frannor gives elements in N(0,1)--ie. pos & negs */
		for (i = 0; i < nfloats; ++i)  noise[i] = frannor();
	break;
	case FLAT: /* franuni gives elements in [0, 1] */
		for (i = 0; i < nfloats; ++i)  noise[i] = 2.0*franuni() - 1.0;
	break;
	default:	/* defensive programming */
		throw cseis_geolib::csException("%d: mysterious itype = %d", __LINE__, itype);
	}


	/* Band limit noise traces if user parObj.getpars any of the f's */
	parObj.getparstring("f",&f);   /* get filter frequencies */
	parObj.getparstring("amps",&amps);  /* get filter amplitudes */
	if ( (*f !='\0') || (*amps !='\0') ) {

		/* Set up call to sufilter */
		char cmdbuf[BUFSIZ];	    /* build sufilter command	*/
		FILE *bandinfp;		    /* fp for input file	*/
		FILE *fp;                   /* fp for pipe to sufilter	*/
		unsigned long nsegy = HDRBYTES + databytes;
		char *segybuf = ealloc1(nsegy, 1);


		/* Trap signals so can remove temp files */
		signal(SIGINT,  (void (*) (int)) closefiles);
		signal(SIGQUIT, (void (*) (int)) closefiles);
		signal(SIGHUP,  (void (*) (int)) closefiles);
		signal(SIGTERM, (void (*) (int)) closefiles);


		/* Prepare temporary files to hold traces */
		bandinfp  = etmpfile();
		bandoutfp = efopen(tmpnam(bandoutfile), "w+");

		/* Paste headers on noise traces and put in tmpfile */
		rewind(headerfp);
		{ register int itr;
		  for (itr = 0; itr < ntr; ++itr) {
			efread(&tr, 1, HDRBYTES, headerfp);
			memcpy((void *) tr.data,
				(const void *)(noise + itr*nt), databytes); 
			fputtr(bandinfp, &tr);
		  }
		}

		/* build cmdbuf ; append sufilter command */
		sprintf(cmdbuf, "sufilter >%s", bandoutfile);

		/* loop through command line args and pass on to cmdbuf */
        	for (--argc, ++argv; argc; --argc, ++argv) {
			/* don't pass sn=, seed=, noise= */
			if ( strncmp(*argv,"sn=",3) &&
			     strncmp(*argv,"seed=",5) &&
			     strncmp(*argv,"noise=",6) ) {
                          strcat(cmdbuf, " ");   /* append a space */
                          strcat(cmdbuf, *argv); /* append the next arg */
			}
        	}

		fp = epopen(cmdbuf, "w");
		rewind (bandinfp);
		{ register int itr;
		  for (itr = 0; itr < ntr; ++itr) {
			efread(segybuf, 1, nsegy, bandinfp);
			efwrite(segybuf, 1, nsegy, fp);
		  }
		}
		efclose(bandinfp);
		epclose(fp);

		/* Load bandlimited traces back into noise vector */
		rewind(bandoutfp);
		{ register int itr;
		  for (itr = 0; itr < ntr; ++itr) {
			fgettr(bandoutfp, &tr);
			memcpy((void *) (noise + itr*nt),
					(const void *) tr.data, databytes); 
		  }
		}
		efclose(bandoutfp);
		eremove(bandoutfile);

	} /* End optional bandlimiting */
		


	/* Compute noise power */
	noipow = 0.0;
	{ register int i;
	  for (i = 0; i < nfloats; ++i) {
		register float noiseval = noise[i];
		noipow += noiseval * noiseval;
	  }
	}


	/* Compute noise scale for desired noise/signal ratio */
	absmaxsig /= sqrt(2.0);  /* make it look like a rmsq value   */
	noipow /= nfloats;	 /* make it the square of rmsq value */
        if( absmaxsig != 0.0 ){
	   noiscale = absmaxsig / (sn * sqrt(noipow));
        }else{
           noiscale = 1.0;
        }


	/* Add scaled noise to trace and output sum */
	rewind(headerfp);
	rewind(tracefp);
	{ register int itr;
	  for (itr = 0; itr < ntr; ++itr) {
		register int trshift = itr*nt;
		register int i;

		efread(&tr, 1, HDRBYTES, headerfp);
		efread(tr.data, 1, databytes, tracefp);
		for (i = 0; i < nt; ++i)
			tr.data[i] += noiscale * noise[trshift + i];

		su2cs->putTrace(&tr);
	  }
	}

	/* Clean up */
	efclose(headerfp);
	if (istmpdir) eremove(headerfile);
	efclose(tracefp);
	if (istmpdir) eremove(tracefile);

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

/* for graceful interrupt termination */
static void closefiles(void)
{
	efclose(headerfp);
	efclose(tracefp);
	eremove(headerfile);
	eremove(tracefile);
	efclose(bandoutfp);
	eremove(bandoutfile);
	exit(EXIT_FAILURE);
}

} // END namespace
