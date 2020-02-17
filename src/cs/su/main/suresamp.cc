/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SURESAMP: $Revision: 1.16 $ ; $Date: 2011/11/16 23:21:55 $        */

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
std::string sdoc_suresamp =
"                                                                   "
" SURESAMP - Resample in time                                       "
"                                                                   "
" suresamp <stdin >stdout  [optional parameters]                    "
"                                                                   "
" Required parameters:                                              "
"     none                                                          "
"                                                                   "
" Optional Parameters:                                              "
"    nt=tr.ns    number of time samples on output                   "
"    dt=         time sampling interval on output                   "
"                default is:                                        "
"                tr.dt/10^6     seismic data                        "
"                tr.d1          non-seismic data                    "
"    tmin=       time of first sample in output                     "
"                default is:                                        "
"                tr.delrt/10^3  seismic data                        "
"                tr.f1          non-seismic data                    "
"    rf=         resampling factor;                                 "
"                if defined, set nt=nt_in*rf and dt=dt_in/rf        "
"    verbose=0   =1 for advisory messages                           "
"                                                                   "
"                                                                   "
" Example 1: (assume original data had dt=.004 nt=256)              "
"    sufilter <data f=40,50 amps=1.,0. |                            "
"    suresamp nt=128 dt=.008 | ...                                  "
" Using the resampling factor rf, this example translates to:       "
"    sufilter <data f=40,50 amps=1.,0. | suresamp rf=0.5 | ...      "
"                                                                   "
" Note the typical anti-alias filtering before sub-sampling!        "
"                                                                   "
" Example 2: (assume original data had dt=.004 nt=256)              "
"    suresamp <data nt=512 dt=.002 | ...                            "
" or use:                                                           "
"    suresamp <data rf=2 | ...                                      "
"                                                                   "
" Example 3: (assume original data had d1=.1524 nt=8192)            "
"    sufilter <data f=0,1,3,3.28 amps=1,1,1,0 |                     "
"    suresamp <data nt=4096 dt=.3048 | ...                          "
"                                                                   "
" Example 4: (assume original data had d1=.5 nt=4096)               "
"    suresamp <data nt=8192 dt=.25 | ...                            "
"                                                                   "
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suresamp {


/* Credits:
 *    CWP: Dave (resamp algorithm), Jack (SU adaptation)
 *    CENPET: Werner M. Heigl - modified for well log support
 *    RISSC: Nils Maercklin 2006 - minor fixes, added rf option
 *
 * Algorithm:
 *    Resampling is done via 8-coefficient sinc-interpolation.
 *    See "$CWPROOT/src/cwp/lib/intsinc8.c" for technical details.
 *
 * Trace header fields accessed:  ns, dt, delrt, d1, f1, trid
 * Trace header fields modified:  ns, dt, delrt (only when set tmin)
 *                                d1, f1 (only when set tmin)
 */
/************************ end self doc ***************************/


segy intrace, outtrace;

void* main_suresamp( void* args )
{
    int nt;            /* number of samples on output trace */
    int nt_in;         /* ... on input trace */
    float dt;          /* sample rate on output trace */
    int idt;           /* ... as integer */
    float dt_in;       /* ... on input trace */
    float tmin;        /* first time sample on output trace */
    float tmin_in;     /* ... on input trace */
    float *t;          /* array of output times */
    int tmin_is_set=0; /* flag for user-defined tmin */
    float rf;          /* resampling factor */
                       /* (rf>1 means down- and rf<1 up-sampling) */ 
    int verbose;       /* if 1(yes) display advisory messages */
    cwp_Bool seismic;  /* flag: is this seismic data? */
 
    
    /* Hook up getpar */
    cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
    cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
    cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
    int argc = suArgs->argc;
    char **argv = suArgs->argv;
    cseis_su::csSUGetPars parObj;

    void* retPtr = NULL;  /*  Dummy pointer for return statement  */
    su2cs->setSUDoc( sdoc_suresamp );
    if( su2cs->isDocRequestOnly() ) return retPtr;
    parObj.initargs(argc, argv);

    try {  /* Try-catch block encompassing the main function body */


    /* Get verbose parameter */
    if (!parObj.getparint("verbose", &verbose)) verbose = 0;

    /* Get information from first trace */
    if (!cs2su->getTrace(&intrace)) throw cseis_geolib::csException("can't get first trace");
    nt_in = intrace.ns;
    if (!nt_in)  throw cseis_geolib::csException("ns not set in header");

    /* check for seismic or well log data */
    seismic = ISSEISMIC(intrace.trid);        
    if (seismic) {
        if (verbose)
            warn("input is seismic data, trid=%d",intrace.trid);
        dt_in   = ((double) intrace.dt)/1000000.0;
        tmin_in = ((double) intrace.delrt)/1000.0;  
    }
    else {
        if (verbose)
            warn("input is not seismic data, trid=%d",intrace.trid);
        dt_in   = intrace.d1;
        tmin_in = intrace.f1;
    }
    
    /* check input times */
    if (!dt_in) parObj.getparfloat("dt_in", &dt_in);
    if (!dt_in) throw cseis_geolib::csException("dt or d1 not set in header or not parObj.getparred");
    if (!tmin_in && verbose) warn("delrt or f1 not set in header");

    /* Get parameters */
    if (!parObj.getparfloat("rf", &rf)) rf=0.0;
    if (rf<0.0) throw cseis_geolib::csException("factor rf=%g must be positive", rf);

    if (rf) {
        if (!parObj.getparint("nt", &nt)) nt = NINT( ((float)nt_in)*rf);
        if (!parObj.getparfloat("dt", &dt)) dt = dt_in/rf;
    }
    else {
        if (!parObj.getparint("nt", &nt)) nt = nt_in;
        if (!parObj.getparfloat("dt", &dt)) dt = dt_in;
    }
    if (parObj.getparfloat("tmin", &tmin)) tmin_is_set = 1;
    parObj.checkpars();

    /* Validate user input nt and dt */
    CHECK_NT("nt",nt);
    idt = NINT(dt * 1000000.0);
        
    /* Allocate vector of output times */
    t = ealloc1float(nt);

    /* Loop over traces */    
    do {
        if (!tmin_is_set)    tmin = tmin_in;
            
        /* Compute output times */
        { register int itime;
          register float tvalue;
          for (itime=0,tvalue=tmin; itime<nt; itime++,tvalue+=dt)
              t[itime] = tvalue;
        }
    
        /* copy and adjust header */
        memcpy(&outtrace, &intrace, HDRBYTES);
        outtrace.ns = nt;
        if (seismic)
            outtrace.dt = idt;
        else
            outtrace.d1 = dt;

        if (tmin_is_set) {
            if (seismic)
                outtrace.delrt = NINT(tmin * 1000.0);
            else
                outtrace.f1 = tmin;
        }
        
        /* sinc interpolate new data */
        ints8r(nt_in, dt_in, tmin_in, intrace.data, 
                0.0, 0.0, nt, t, outtrace.data);
        
        su2cs->putTrace(&outtrace);
    } while (cs2su->getTrace(&intrace));


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
