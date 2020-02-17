/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUTAUP: $Revision: 1.14 $ ; $Date: 2011/11/16 23:35:04 $	*/

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
#include "taup.h"
extern "C" {
  #include <signal.h>

}

/*********************** self documentation **********************/
std::string sdoc_sutaup =
"                                                                       "
" SUTAUP - forward and inverse T-X and F-K global slant stacks		"
"                                                                       "
"    sutaup <infile >outfile  [optional parameters]                 	"
"                                                                       "
" Optional Parameters:                                                  "
" option=1			=1 for forward F-K domian computation	"
"				=2 for forward T-X domain computation	"
"				=3 for inverse F-K domain computation	"
"				=4 for inverse T-X domain computation	"
" dt=tr.dt (from header) 	time sampling interval (secs)           "
" nx=ntr   (counted from data)	number of horizontal samples (traces)	"
" dx=1				horizontal sampling interval (m)	"
" npoints=71			number of points for rho filter		"
" pmin=0.0			minimum slope for Tau-P transform (s/m)	"
" pmax=.006			maximum slope for Tau-P transform (s/m)	"
" np=nx				number of slopes for Tau-P transform	"
" ntau=nt			number of time samples in Tau-P domain  "
" fmin=3			minimum frequency of interest 	        "
" xmin=0			offset on first trace	 	        "
"                                                                       "
" verbose=0	verbose = 1 echoes information				"
"									"
" tmpdir= 	 if non-empty, use the value as a directory path	"
"		 prefix for storing temporary files; else if the	"
"	         the CWP_TMPDIR environment variable is set use		"
"	         its value for the path; else use tmpfile()		"
" 									"
" Notes:                                                                "
" The cascade of a forward and inverse  tau-p transform preserves the	"
" relative amplitudes in a data panel, but not the absolute amplitudes  "
" meaning that a scale factor must be applied to data output by such a  "
" a cascade before the output may be compared to the original data.	"
" This is a characteristic of the algorithm employed in this program.	"
" (Suradon does not have this problem.)					"
"                                                                       "
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sutaup {


/*
 * Credits: CWP: Gabriel Alvarez, 1995.
 *
 * Reference:       
 *    Levin, F., editor, 1991, Slant-Stack Processing, Geophysics Reprint 
 *         Series #14, SEG Press, Tulsa.
 *
 * Trace header fields accessed: ns, dt
 * Trace header fields modified: dt,d2,f2
 */
/**************** end self doc ********************************/

static void closefiles(void);

/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/


segy tr;

void* main_sutaup( void* args )
{
	int ix,it;		/* loop counters */
	int ntr;		/* number of input traces */
	int nt;			/* number of time samples */
	int nx;			/* number of horizontal samples */
	int option;		/* flag for requested opeartion */
	float dt;               /* Time sample interval */
        float dx;               /* horizontal sample interval */
	float xmin;		/* offset on first trace */
        float pmin;             /* Minimum slope for Tau-P transform */
        float pmax;             /* Maximum slope for Tau-P transform */
	float dp;		/* slope sampling interval */
	int np;			/* number of slopes for slant stack */
	float fmin;		/* minimum frequency of interest */
	int npoints;		/* number of points for rho filter */
	float **in_traces=NULL;	/* array[nx][nt] of input traces */	
	float **out_traces=NULL;/* array[nx][nt] of output traces */	
	int verbose;		/* flag for echoing information */
	char *tmpdir=NULL;	/* directory path for tmp files */
	cwp_Bool istmpdir=cwp_false;/* true for user-given path */
	
        /* hook up getpar to handle the parameters */
        cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
        cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
        cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
        int argc = suArgs->argc;
        char **argv = suArgs->argv;
        cseis_su::csSUGetPars parObj;

        void* retPtr = NULL;  /*  Dummy pointer for return statement  */
        su2cs->setSUDoc( sdoc_sutaup );
        if( su2cs->isDocRequestOnly() ) return retPtr;
        parObj.initargs(argc, argv);

        try {  /* Try-catch block encompassing the main function body */


	if (!parObj.getparint("verbose", &verbose))	verbose = 0;

	/* Look for user-supplied tmpdir */
	if (!parObj.getparstring("tmpdir",&tmpdir) &&
	    !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
	if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
		throw cseis_geolib::csException("you can't write in %s (or it doesn't exist)", tmpdir);

        /* get info from first trace */
        if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
        nt = tr.ns;
        dt = (float) tr.dt/1000000.0;

        /* Store traces in tmpfile while getting a count */
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
        ntr = 0;
        do {
                ++ntr;
                efwrite(&tr, 1, HDRBYTES, headerfp);
                efwrite(tr.data, FSIZE, nt, tracefp);
        } while (cs2su->getTrace(&tr));

        /* get general flags and parameters and set defaults */
        if (!parObj.getparint("np",&np))             	np = ntr;
        if (!parObj.getparint("nx",&nx))          	nx = ntr;
        if (!parObj.getparint("option",&option))       option = 1;
        if (!parObj.getparfloat("pmin",&pmin))		pmin = 0.0;
        if (!parObj.getparfloat("pmax",&pmax))		pmax = 0.006;
        if (!parObj.getparfloat("xmin",&xmin))		xmin = 0.0;
	if (!parObj.getparfloat("dx",&dx))		dx = 1.0;
	if (!parObj.getparfloat("dt",&dt))		dt = dt;
	if (!parObj.getparfloat("fmin",&fmin))		fmin = 3.;
	if (!parObj.getparint("npoints",&npoints))	npoints = 71;


        parObj.checkpars();

	if (dt == 0.0)
		throw cseis_geolib::csException("header field dt not set, must be parObj.getparred");

	/* allocate space */
        in_traces = alloc2float(nt, ntr);
        out_traces = alloc2float(nt, np);
	dp=(pmax-pmin)/(np-1);

        /* load traces into an array and close temp file */
	erewind(headerfp);
        erewind(tracefp);
        for (ix=0; ix<ntr; ix++)
                fread (in_traces[ix], FSIZE, nt, tracefp);
        efclose (tracefp);
	if (istmpdir) eremove(tracefile);

	/* do requested operation */ 
	if (option==1) {
		/* compute F-K forward slant stack */
		fwd_FK_sstack (dt, nt, nx, xmin, dx, np, pmin, dp, fmin,
	        	in_traces, out_traces);
	} else if (option==2) {
		/* compute t-x forward slant stack */
		fwd_tx_sstack (dt, nt, nx, xmin, dx, np, pmin, dp, 
	        	in_traces, out_traces);
	} else if (option==3) {
		/* compute F-K inverse slant stack */
		inv_FK_sstack (dt, nt, nx, xmin, dx, np, pmin, dp, fmin,
	        	in_traces, out_traces);
	} else if (option==4) {
		/* compute t-x inverse slant stack */
		inv_tx_sstack (dt, nt, nx, npoints, xmin, dx, np, pmin, dp,
			in_traces, out_traces);
	} else throw cseis_geolib::csException("option flag has to be between 1 and 4");
		
        /* write output traces */
        erewind(headerfp);
	{       register int itr;
		for (itr=0; itr<np; itr++) {
			efread(&tr, 1, HDRBYTES, headerfp);
                        tr.tracl = 1+itr;
                        tr.tracr = 1+itr;
                        tr.dt=(int)(dt*1000000.0);
                        tr.ns=nt;
                        tr.d2 = dp;
                        tr.f2 = pmin;
            
			for (it=0; it<nt; it++) 
				tr.data[it]=out_traces[itr][it];
			su2cs->putTrace(&tr);
		}
	}
	efclose(headerfp);
	if (istmpdir) eremove(headerfile);

	/* free allocated space */
	free2float(in_traces);
	free2float(out_traces);

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
	exit(EXIT_FAILURE);
}

} // END namespace
