/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUTVBAND: $Revision: 1.18 $ ; $Date: 2011/11/12 00:09:00 $      */

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
std::string sdoc_sutvband =
"                                                               "
" SUTVBAND - time-variant bandpass filter (sine-squared taper)  "
"                                                               "
" sutvband <stdin >stdout tf= f=			        "
"                                                               "
" Required parameters:                                          "
"       dt = (from header)      time sampling interval (sec)    "
"       tf=             times for which f-vector is specified   "
"       f=f1,f2,f3,f4   Corner frequencies corresponding to the "
"                       times in tf. Specify as many f= as      "
"                       there are entries in tf.                "
"                                                               "
" The filters are applied in frequency domain.                  "
"                                                               "
" Example:                                                      "
" sutvband <data tf=.2,1.5 f=10,12.5,40,50 f=10,12.5,30,40 | ..."
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sutvband {


/* Credits:
 *      CWP: Jack, Ken
 *
 * Trace header fields accessed:  ns, dt, delrt
 */
/**************** end self doc ***********************************/


#define PIBY2   1.57079632679490
#define LOOKFAC 2       /* Look ahead factor for npfaro   */
#define PFA_MAX 720720  /* Largest allowed nfft           */

/* Prototypes */
void makefilter(float *f, int nfft, int nfreq, float dt, float *filter);
void bandpass(float *data, int ntime, int nfft, int nfreq, 
		float *filterj, float *ftracej);

segy tr;

void* main_sutvband( void* args )
{
        float **filter;         /* filter arrays                        */
        float *tf=NULL;		/* times at which filters are centered  */
        int *itf;               /* ... as integers                      */
	int jmin;		/* index of first filter itf value	*/
	int jmax;		/* index of last filter itf value	*/
        int nfft;     	        /* fft sizes in each time gate          */
        int nfreq;     	        /* number of frequencies  	        */
        float **ftrace;         /* filtered sub-traces                  */
        int nfilter;            /* number of filters specified          */
      	float dt;               /* sample spacing                       */
        float tmin;             /* first time on traces                 */
        int nt;                 /* number of points on input trace      */
	cwp_Bool seismic;	/* is this seismic data?		*/
	float *data;

        
        /* Initialize */
        cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
        cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
        cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
        int argc = suArgs->argc;
        char **argv = suArgs->argv;
        cseis_su::csSUGetPars parObj;

        void* retPtr = NULL;  /*  Dummy pointer for return statement  */
        su2cs->setSUDoc( sdoc_sutvband );
        if( su2cs->isDocRequestOnly() ) return retPtr;
        parObj.initargs(argc, argv);

        try {  /* Try-catch block encompassing the main function body */



        /* Get info from first trace */ 
        if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	seismic = ISSEISMIC(tr.trid); 
		 
	if (!seismic)
		warn("input is not seismic data, trid=%d", tr.trid);
        nt = tr.ns;
        if (!parObj.getparfloat("dt", &dt))    dt = ((double) tr.dt)/1000000.0;
        if (!dt) throw cseis_geolib::csException("dt field is zero and not parObj.getparred");
	tmin = tr.delrt/1000.0;


        /* Get number of filters and center times */
        if (!(nfilter = parObj.countparval("tf")))  CSMUSTGETPARFLOAT("tf", tf);
	if (parObj.countparname("f") != nfilter)
		throw cseis_geolib::csException("must give one f 4-tuple for each"
		    " (%d) tf value", nfilter);
		
	/* Leave room for possibly missing filters at endpoints */
	tf = ealloc1float(nfilter+4);  /* never use first 2 or last 2 */
	itf = ealloc1int(nfilter+4);
        parObj.getparfloat("tf", tf+2); jmin = 2; jmax = nfilter + 1;
	{ register int j;
          for (j = jmin; j <= jmax; ++j)  itf[j] = NINT((tf[j] - tmin)/dt);
        }
	

        /* Make filters with scale for inverse transform */
	nfft = npfaro(nt, LOOKFAC * nt);
	if (nfft >= SU_NFLTS || nfft >= PFA_MAX)
		throw cseis_geolib::csException("Padded nt=%d -- too big", nfft);
	nfreq = nfft/2 + 1;
	filter = ealloc2float(nfreq, nfilter+4); /* never use 1st & last */
        { register int j;
          for (j = jmin; j <= jmax; ++j) {
                float *f = ealloc1float(4);

                if (parObj.getnparfloat(j-jmin+1, "f", f) != 4)
                        throw cseis_geolib::csException("must give 4 corner frequencies in f=");
        	if (f[0] < 0.0 || f[0] > f[1] ||
				  f[1] >= f[2] || f[2] > f[3])
               		throw cseis_geolib::csException("Filter #%d has bad frequencies", j - jmin + 1);
                makefilter(f, nfft, nfreq, dt, filter[j]);
          }
        }
	

	/* User may not have given a filter for tmin and/or tmax--	*/
	/* Extend array so can always assume these filters are present.	*/
	/* Note don't really use any of the extra storage in **filter!	*/
	if (itf[jmin] > 0) {
		filter[jmin-1] = filter[jmin]; 
		itf[jmin-1] = 0;
		--jmin;
	}
	if (itf[jmax] < nt - 1) {
		filter[jmax+1] = filter[jmax];
		itf[jmax+1] = nt - 1;
		++jmax;
	}
	
	
	/* Extend array so can always consider time points to be interior */
	itf[jmin-1] = 0;      /* now jmin - 1 is a valid index */
	itf[jmax+1] = nt - 1; /* now jmax + 1 is a valid index */


        /* Main loop over traces */
        ftrace = ealloc2float(nt, nfilter+4); /* never use 1st & last */
	data = ealloc1float(nt);
        do {
                register int i, j;
		
		/* Construct filtered sub-traces */
		for (j = jmin; j <= jmax; ++j) {			
			memset( (void *) data, 0, nt*FSIZE);
			for (i = itf[j-1]; i <= itf[j+1]; ++i)
				data[i] = tr.data[i];
                        bandpass(data,nt,nfft,nfreq,filter[j],ftrace[j]);
                }
		               
               /* Compose filtered trace from sub-traces */
               for (j = jmin; j < jmax; ++j) {
	       		float fitfj;
                        for (fitfj = i = itf[j]; i <= itf[j+1]; ++i) {
                                float a = (i - fitfj)/(itf[j+1] - fitfj);
                                tr.data[i] = (1-a)*ftrace[j][i] +
                                                 a*ftrace[j+1][i];
			}
                }
                
                su2cs->putTrace(&tr);
        } while (cs2su->getTrace(&tr));

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

void bandpass(float *data, int nt, int nfft, int nfreq,
		float *filterj, float *ftracej)
{
	float *rt = ealloc1float(nfft);
	complex *ct = ealloc1complex(nfreq);
	register int i;

        /* Load trace into rt (zero-padded) */
        memcpy( (void *) rt, (const void *) data, nt*FSIZE);
        memset( (void *) (rt + nt), 0, (nfft-nt)*FSIZE);

        /* FFT, filter, inverse FFT */
        pfarc(1, nfft, rt, ct);
        for (i = 0; i < nfreq; ++i)  ct[i] = crmul(ct[i], filterj[i]);
        pfacr(-1, nfft, ct, rt);

        /* Load traces back in, recall filter had nfft factor */
        for (i = 0; i < nt; ++i)  ftracej[i] = rt[i]; /* ftracej = rt ?? */

	/* free temp arrays */
	free1float( rt );
	free1complex( ct );
}


void makefilter(float *f, int nfft, int nfreq, float dt, float *filter)
{
        float onfft = 1.0 / nfft;
        float df = onfft / dt;
        int nfreqm1 = nfreq - 1;
        int if1 = NINT(f[0]/df);
        int if2 = NINT(f[1]/df);
        int if3 = MIN(NINT(f[2]/df), nfreqm1);
        int if4 = MIN(NINT(f[3]/df), nfreqm1);


        { register int i;
	  register float c = PIBY2 / (if2 - if1 + 2);
	  for (i = if1; i <= if2; ++i) {
		register float s = sin(c*(i - if1 + 1));
		filter[i] = s * s * onfft;
          }
	 }

        { register int i;
	  register float c = PIBY2 / (if4 - if3 + 2);
	  for (i = if3; i <= if4; ++i) {
		register float s = sin(c*(if4 - i + 1));
		filter[i] = s * s * onfft;
	  }
        }

        { register int i;
          for (i = if2 + 1; i < if3;   ++i)  filter[i] = onfft; 
          for (i = 0;       i < if1;   ++i)  filter[i] = 0.0; 
          for (i = if4 + 1; i < nfreq; ++i)  filter[i] = 0.0; 
        }
}

} // END namespace
