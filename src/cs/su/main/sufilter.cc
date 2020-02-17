/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUFILTER: $Revision: 1.23 $ ; $Date: 2011/11/12 00:09:00 $        */


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
std::string sdoc_sufilter =
"									"
" SUFILTER - applies a zero-phase, sine-squared tapered filter		"
"									"
" sufilter <stdin >stdout [optional parameters]         		"
"									"
" Required parameters:                                         		"
"       if dt is not set in header, then dt is mandatory        	"
"									"
" Optional parameters:							"
"       f=f1,f2,...             array of filter frequencies(HZ) 	"
"       amps=a1,a2,...          array of filter amplitudes		"
"       dt = (from header)      time sampling interval (sec)        	"
"	verbose=0		=1 for advisory messages		"
"									"
" Defaults:f=.10*(nyquist),.15*(nyquist),.45*(nyquist),.50*(nyquist)	"
"                        (nyquist calculated internally)		"
"          amps=0.,1.,...,1.,0.  trapezoid-like bandpass filter		"
"									"
" Examples of filters:							"
" Bandpass:   sufilter <data f=10,20,40,50 | ...			"
" Bandreject: sufilter <data f=10,20,30,40 amps=1.,0.,0.,1. | ..	"
" Lowpass:    sufilter <data f=10,20,40,50 amps=1.,1.,0.,0. | ...	"
" Highpass:   sufilter <data f=10,20,40,50 amps=0.,0.,1.,1. | ...	"
" Notch:      sufilter <data f=10,12.5,35,50,60 amps=1.,.5,0.,.5,1. |..	"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sufilter {


/* Credits:
 *      CWP: John Stockwell, Jack Cohen
 *	CENPET: Werner M. Heigl - added well log support
 *
 * Possible optimization: Do assignments instead of crmuls where
 * filter is 0.0.
 *
 * Trace header fields accessed: ns, dt, d1
 */
/**************** end self doc ***********************************/

/* Prototype of function used internally */
void polygonalFilter(float *f, float *amps,
			int npoly, int nfft, float dt, float *filter);

#define PIBY2   1.57079632679490
#define FRAC0   0.10    /* Ratio of default f1 to Nyquist */
#define FRAC1   0.15    /* Ratio of default f2 to Nyquist */
#define FRAC2   0.45    /* Ratio of default f3 to Nyquist */
#define FRAC3   0.50    /* Ratio of default f4 to Nyquist */
#define LOOKFAC 2       /* Look ahead factor for npfao    */
#define PFA_MAX 720720  /* Largest allowed nfft           */


segy tr;

void* main_sufilter( void* args )
{
        register float *rt;     /* real trace                           */
        register complex *ct;   /* complex transformed trace            */
        float *filter;          /* filter array                         */
        float *f;               /* array of filter frequencies		*/
        int npoly;              /* .... sizes of f and intfr	        */
        float *amps;            /* array of amplitude values		*/
        int namps;              /* .... size of amps                    */
        int icount,ifs,iamps;   /* loop counting variables              */
        float dt;               /* sample spacing                       */
        float nyq;              /* nyquist frequency                    */
        int nt;                 /* number of points on input trace      */
        int nfft;               /* number of points for fft trace       */
        int nf;                 /* number of frequencies (incl Nyq)     */
	int verbose;		/* flag to get advisory messages	*/
	cwp_Bool seismic;	/* is this seismic data?		*/

        
        /* Initialize */
        cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
        cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
        cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
        int argc = suArgs->argc;
        char **argv = suArgs->argv;
        cseis_su::csSUGetPars parObj;

        void* retPtr = NULL;  /*  Dummy pointer for return statement  */
        su2cs->setSUDoc( sdoc_sufilter );
        if( su2cs->isDocRequestOnly() ) return retPtr;
        parObj.initargs(argc, argv);

        try {  /* Try-catch block encompassing the main function body */



        /* Get info from first trace */ 
	if (!parObj.getparint("verbose", &verbose))	verbose=0;
        if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	seismic = ISSEISMIC(tr.trid);
	if (seismic) {
		if (verbose)	warn("input is seismic data, trid=%d",tr.trid);
		dt = ((double) tr.dt)/1000000.0;
	}
	else {
		if (verbose) warn("input is not seismic data, trid=%d",tr.trid);

		dt = tr.d1;

        }

	/* error trapping so that the user can have a default value of dt */
	if (!(dt || parObj.getparfloat("dt", &dt))) {
		dt = .004;
		if (verbose) {
			warn("neither dt nor d1 are set, nor is dt parObj.getparred!");
			warn("assuming .004 sec sampling!");
		}
	}

        nt = tr.ns;
        nyq = 0.5/dt;


        /* Set up FFT parameters */
        nfft = npfaro(nt, LOOKFAC * nt);
        if (nfft >= SU_NFLTS || nfft >= PFA_MAX)
                throw cseis_geolib::csException("Padded nt=%d -- too big", nfft);

        nf = nfft/2 + 1;

        /* Get frequencies that define the filter */
        if ((npoly = parObj.countparval("f"))!=0) {
                f = ealloc1float(npoly);
                parObj.getparfloat("f",f);
        } else {
                npoly = 4;
                f = ealloc1float(npoly);

                f[0] = FRAC0 * nyq;
                f[1] = FRAC1 * nyq;
                f[2] = FRAC2 * nyq;
                f[3] = FRAC3 * nyq;
        }

	/* Check f values */
	if(npoly < 2) warn("Only %d value defining filter",npoly);
        for(ifs=0; ifs < npoly-1; ++ifs)
		if(f[ifs] < 0.0 || f[ifs] > f[ifs+1])
                                throw cseis_geolib::csException("Bad filter parameters");
	
	/* Get filter amplitude values*/
        if ((namps = parObj.countparval("amps"))!=0) {
                amps = ealloc1float(namps);
                parObj.getparfloat("amps",amps);
        } else {
                namps = npoly;
                amps = ealloc1float(namps);

		/* default is a trapezoidal bandpass filter */
		for(iamps=0; iamps<namps; ++iamps)
               		amps[iamps]=1.;
		
		amps[0]=0.; amps[namps-1]=0.;
        }
	if (!(namps==npoly)) 
		throw cseis_geolib::csException("number of f values must = number of amps values");
        
        /* Check amps values */
        for(iamps = 0, icount=0; iamps < namps ; ++iamps) {
		if( amps[iamps] > 0. ) ++icount;
                if( amps[iamps] < 0.) throw cseis_geolib::csException("amp values must be positive");
        }
        if (icount==0) throw cseis_geolib::csException("All amps values are zero");
        for(iamps = 0, icount=0; iamps < namps-1 ; ++iamps) {
			if(!(amps[iamps]==amps[iamps+1])) ++icount;
	}
        if (icount==0) warn("All amps values are the same");


        /* Allocate fft arrays */
        rt   = ealloc1float(nfft);
        ct   = ealloc1complex(nf);
        filter = ealloc1float(nf);

	/* Build the polygonal filter filter[]*/
	polygonalFilter(f,amps,npoly,nfft,dt,filter);

        /* Main loop over traces */
        do {
                register int i;

                /* Load trace into rt (zero-padded) */
                memcpy((void *) rt, (const void *) tr.data, nt*FSIZE);
                memset((void *) (rt + nt), 0 , (nfft-nt)*FSIZE);

                /* FFT, filter, inverse FFT */
                pfarc(1, nfft, rt, ct);
                for (i = 0; i < nf; ++i)  ct[i] = crmul(ct[i], filter[i]);
                pfacr(-1, nfft, ct, rt);

                /* Load traces back in, recall filter had nfft factor */
                for (i = 0; i < nt; ++i)  tr.data[i] = rt[i];

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


void polygonalFilter(float *f, float *amps, int npoly,
				int nfft, float dt, float *filter)
/*************************************************************************
polygonalFilter -- polygonal filter with sin^2 tapering
**************************************************************************
Input:
f		array[npoly] of frequencies defining the filter
amps		array[npoly] of amplitude values
npoly		size of input f and amps arrays
dt		time sampling interval
nfft		number of points in the fft

Output:
filter		array[nfft] filter values
**************************************************************************
Notes: Filter is to be applied in the frequency domain
**************************************************************************
Author:  CWP: John Stockwell   1992
*************************************************************************/
#define PIBY2   1.57079632679490
{
        int *intfr;             /* .... integerizations of f		*/
        int icount,ifs;		/* loop counting variables              */
	int taper=0;		/* flag counter				*/
        int nf;                 /* number of frequencies (incl Nyq)     */
        int nfm1;               /* nf-1                                 */
        float onfft;            /* reciprocal of nfft                   */
        float df;               /* frequency spacing (from dt)          */

        
	intfr=alloc1int(npoly);

        nf = nfft/2 + 1;
        nfm1 = nf - 1;
        onfft = 1.0 / nfft;

        /* Compute array of integerized frequencies that define the filter*/
        df = onfft / dt;
        for(ifs=0; ifs < npoly ; ++ifs) {
                intfr[ifs] = NINT(f[ifs]/df);
                if (intfr[ifs] > nfm1) intfr[ifs] = nfm1;
        }

	/* Build filter, with scale, and taper specified by amps[] values*/
	/* Do low frequency end first*/
	for(icount=0; icount < intfr[0] ; ++icount) 
		filter[icount] = amps[0] * onfft;

	/* now do the middle frequencies */
	for(ifs=0 ; ifs<npoly-1 ; ++ifs){
	   if(amps[ifs] < amps[ifs+1]) {	
		++taper;
		for(icount=intfr[ifs]; icount<=intfr[ifs+1]; ++icount) {
		    float c = PIBY2 / (intfr[ifs+1] - intfr[ifs] + 2);
		    float s = sin(c*(icount - intfr[ifs] + 1));
		    float adiff = amps[ifs+1] - amps[ifs];
		    filter[icount] = (amps[ifs] + adiff*s*s) * onfft;
		}
	   } else if (amps[ifs] > amps[ifs+1]) {	
		++taper;
		for(icount=intfr[ifs]; icount<=intfr[ifs+1]; ++icount) {
			   float c = PIBY2 / (intfr[ifs+1] - intfr[ifs] + 2);
                	   float s = sin(c*(intfr[ifs+1] - icount + 1));
			   float adiff = amps[ifs] - amps[ifs+1];
                	   filter[icount] = (amps[ifs+1] + adiff*s*s) * onfft;
		  }
	   } else 
		if(!(taper)){
		for(icount=intfr[ifs]; icount <= intfr[ifs+1]; ++icount)
		   	   filter[icount] = amps[ifs] * onfft;
		} else {
		for(icount=intfr[ifs]+1; icount <= intfr[ifs+1]; ++icount)
		   	   filter[icount] = amps[ifs] * onfft;
		}
	}

	/* finally do the high frequency end */
	for(icount=intfr[npoly-1]+1; icount<nf; ++icount){
		filter[icount] = amps[npoly-1] * onfft;
	}

}

} // END namespace
