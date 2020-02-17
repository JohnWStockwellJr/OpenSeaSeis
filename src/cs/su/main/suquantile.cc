/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUQUANTILE: $Revision: 1.9 $ ; $Date: 2011/11/16 17:24:58 $         */

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
extern "C" {
  #include <signal.h>

}

/*********************** self documentation *****************************/
std::string sdoc_suquantile =
"                                                                       "
" SUQUANTILE - display some quantiles or ranks of a data set            "
"                                                                       "
" suquantile <stdin >stdout [optional parameters]			"
"                                                                       "
" Required parameters:                                                  "
"       none (no-op)                                                    "
"                                                                       "
" Optional parameters:                                                  "
"	panel=1		flag; 0 = do trace by trace (vs. whole data set)"
"	quantiles=1	flag; 0 = give ranks instead of quantiles	"
" 	verbose=0	verbose = 1 echoes information			"
"									"
" 	tmpdir= 	 if non-empty, use the value as a directory path"
"			 prefix for storing temporary files; else if the"
"		         the CWP_TMPDIR environment variable is set use	"
"		         its value for the path; else use tmpfile()	"
" 									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suquantile {


/* Credits:
 *      CWP: Jack K. Cohen
 *
 *
 * Trace header fields accessed: ns, tracl, mark
 */
/**************** end self doc *******************************************/

/* subroutine prototypes */
static void get_quantiles(int *rank, int nt, int fsize);
static void get_ranks(int *rank, int nt, int fsize);
int cmp_indirect(int *r, int *s);
static void closefiles(void);

/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
FILE *tracefp;		/* fp for trace storage file		*/

float *data;		/* the data; global to use system qsort */

segy tr;

void* main_suquantile( void* args )
{
	int quantiles;	/* flag for rank instead of quantiles	*/
	int panel;	/* flag to report trace by trace	*/
	int nt;		/* number of samples on trace		*/
	int *rank;	/* permuted indices for indirect sort	*/
	int verbose;	/* flag for echoing information		*/
	char *tmpdir;	/* directory path for tmp files		*/
	cwp_Bool istmpdir=cwp_false;/* true for user-given path	*/


        /* Initialize */
        cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
        cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
        cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
        int argc = suArgs->argc;
        char **argv = suArgs->argv;
        cseis_su::csSUGetPars parObj;

        void* retPtr = NULL;  /*  Dummy pointer for return statement  */
        su2cs->setSUDoc( sdoc_suquantile );
        if( su2cs->isDocRequestOnly() ) return retPtr;
        parObj.initargs(argc, argv);

        try {  /* Try-catch block encompassing the main function body */


        /* Get nt from first trace */
        if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
        nt   = (int) tr.ns;

        /* Get parameters */
	if (!parObj.getparint("panel", &panel))  panel = 1;
        if (!parObj.getparint("quantiles", &quantiles))  quantiles = 1;
	if (!parObj.getparint("verbose", &verbose))	verbose = 0;

	/* Look for user-supplied tmpdir */
	if (!parObj.getparstring("tmpdir",&tmpdir) &&
	    !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
	if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
		throw cseis_geolib::csException("you can't write in %s (or it doesn't exist)", tmpdir);
        parObj.checkpars();

	/* Main loop over traces */
	if (!panel) { /* trace by trace */
		rank = ealloc1int(nt);	/* array for indirect sort */
		data = ealloc1float(nt);
		fprintf(stderr, "trace length: %d\n", nt);
		do {
			fprintf(stderr, "trace: %d\n", tr.tracl);	
			memcpy((void *)data, (const void *) tr.data, nt*FSIZE);

			if (!tr.mark) {
				if (!quantiles) { /* give ranks */
					get_ranks(rank, nt, FSIZE);
				} else { /* give quantiles */
					get_quantiles(rank, nt, FSIZE);
				}
				memcpy((void *)tr.data, (const void *) data,
				       nt*FSIZE);
			}
		} while(cs2su->getTrace(&tr));
	} else { /* do whole data set at once */
		int ndata, ntr = 0;

		if (STREQ(tmpdir,"")) {
			tracefp = etmpfile();
			if (verbose) warn("using tmpfile() call");
		} else { /* user-supplied tmpdir */
			char directory[BUFSIZ];
			strcpy(directory, tmpdir);
			strcpy(tracefile, temporary_filename(directory));
			/* Trap signals so can remove temp files */
			signal(SIGINT,  (void (*) (int)) closefiles);
			signal(SIGQUIT, (void (*) (int)) closefiles);
			signal(SIGHUP,  (void (*) (int)) closefiles);
			signal(SIGTERM, (void (*) (int)) closefiles);
			tracefp = efopen(tracefile, "w+");
			istmpdir=cwp_true;		
			if (verbose)
			     warn("putting temporary files in %s", directory);
		}

		do {
			++ntr;
			efwrite(tr.data, FSIZE, nt, tracefp);
		} while (cs2su->getTrace(&tr));
		erewind(tracefp);

		ndata = nt*ntr;
		fprintf(stderr, "trace length: %d\n", nt);
		fprintf(stderr, "number of traces: %d\n", ntr);
		fprintf(stderr, "number of samples: %d\n", ndata);
		rank = ealloc1int(ndata);  /* array for indirect sort */
		data = ealloc1float(ndata);
		
		/* Load traces into data and close tmpfile */
		efread(data, FSIZE, ndata, tracefp);
		efclose(tracefp);
		
		if (!quantiles) { /* give ranks */
			get_ranks(rank, ndata, FSIZE);
		} else { /* give quantiles */
			get_quantiles(rank, ndata, FSIZE);
		}
	}
	
	/* Clean up */
	if (istmpdir) eremove(tracefile);
	free1(data);
	
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

/* GET_RANKS - fprint some sample ranks */
static void get_ranks(int *rank, int nt, int fsize)
{
	register int i;
	
	for (i = 0; i < nt; ++i)  rank[i] = i;
	qsort(rank, nt, fsize, (int (*)()) cmp_indirect);
	
	/* Could make table of desired i's and loop */
	i = 0;
	fprintf(stderr, " rank[%d] = %8.2e", i+1, data[rank[i]]);
	i = nt / 20;
	fprintf(stderr, " rank[%d] = %8.2e", i+1, data[rank[i]]);
	i = nt/2 - i;
	fprintf(stderr, " rank[%d] = %8.2e", i+1, data[rank[i]]);
	fprintf(stderr, "\n");
	i = nt - 1 - i;
	fprintf(stderr, " rank[%d] = %8.2e", i+1, data[rank[i]]);
	i = nt - 1 - nt/20;
	fprintf(stderr, " rank[%d] = %8.2e", i+1, data[rank[i]]);
	i = nt - 1;
	fprintf(stderr, " rank[%d] = %8.2e", i+1, data[rank[i]]);
	fprintf(stderr, "\nmin is at sample %d,  max at %d\n",
	       rank[0] + 1, rank[nt-1] + 1);
}

/* GET_QUANTILES - fprint some sample quantiles */
static void get_quantiles(int *rank, int nt, int fsize)
{
	register int i, iq;
	for (i = 0; i < nt; ++i)  rank[i] = i;
	qsort(rank, nt, fsize, (int (*)()) cmp_indirect);

	/* round to qth quantile (100 qth percentile) */
	/* thus (q*nt - 1) + .5 (-1 for zero basing) */
	i = 1; iq = (int) (0.01*nt - 0.5);
	fprintf(stderr, " %dst percentile is %8.2e\n",
	       i+1, data[rank[iq]]);
	i = 5; iq = (int) (0.05*nt - 0.5);
	fprintf(stderr, " %dth percentile is %8.2e\n",
	       i+1, data[rank[iq]]);
	i = 25; iq = (int) (0.25*nt - 0.5);
	fprintf(stderr, "%dth percentile is %8.2e\n",
	       i+1, data[rank[iq]]);
	i = 50; iq = (int) (0.50*nt - 0.5);
	fprintf(stderr, "%dth percentile is %8.2e\n",
	       i+1, data[rank[iq]]);
	i = 75; iq = (int) (0.75*nt - 0.5);
	fprintf(stderr, "%dth percentile is %8.2e\n",
	       i+1, data[rank[iq]]);
	i = 95; iq = (int) (0.95*nt - 0.5);
	fprintf(stderr, "%dth percentile is %8.2e\n",
	       i+1, data[rank[iq]]);
	i = 99; iq = (int) (0.99*nt - 0.5);
	fprintf(stderr, "%dth percentile is %8.2e\n",
	       i+1, data[rank[iq]]);
	fprintf(stderr, "min at sample %d equals %8.2e\n",
	       rank[0] + 1, data[rank[0]]);
	fprintf(stderr, "max at sample %d equals %8.2e\n",
	       rank[nt-1] + 1, data[rank[nt-1]]);
}


/* Comparison function for qsort */
int cmp_indirect(int *r, int *s)
{
	float diff = data[*r] - data[*s];

	if      (diff > 0)	return(1);
	else if (diff < 0)	return(-1);
	else  /* diff == 0 */	return(0);
}

/* for graceful interrupt termination */
static void closefiles(void)
{
	efclose(tracefp);
	eremove(tracefile);
	exit(EXIT_FAILURE);
}

} // END namespace
