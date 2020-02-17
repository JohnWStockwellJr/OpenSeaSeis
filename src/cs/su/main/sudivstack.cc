/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUDIVSTACK: $Revision: 1.12 $ ; $Date: 2011/11/16 23:14:54 $	*/

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
std::string sdoc_sudivstack =
"                                                                       "
" SUDIVSTACK -  Diversity Stacking using either average power or peak   "
"               power within windows                                    "
"                                                                       "
"                                                                       "
" Required parameters:                                                  "
"    none                                                               "
"                                                                       "
" Optional parameters:                                                  "
"    key=tracf        key header word to stack on                       "
"    winlen=0.064     window length in seconds.                         "
"                     typical choices: 0.064, 0.128, 0.256,             "
"                                      0.512, 1.024, 2.048, 4.096       "
"                     if not specified the entire trace is used         "
"                                                                       "
"    peak=1           peak power option default is average power        "
"                                                                       "
" Notes:                                                                "
"    Diversity stacking is a noise reduction technique used in the      "
"    summation of duplicate data. Each trace is scaled by the inverse   "
"    of its average power prior to stacking.  The composite trace is    "
"    then renormalized by dividing by the sum of the scalers used.      "
"                                                                       "
"    This program stacks adjacent traces having the same key header     "
"    word, which can be specified by the key parameter. The default     "
"    is \"tracf\" (trace number within field record).                   "
"    For more information on key header words, type \"sukeyword -o\".   "
"                                                                       "
" Examples:                                                             "
"    For duplicate field shot records:                                  "
"        susort < field.data tracf | sudivstack > stack.data            "
"    For CDP ordered data:                                              "
"        sudivstack < cdp.data key=cdp > stack.data                     "
"                                                                       "
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sudivstack {


/* 
 * Author: Mary Palen-Murphy,
 *         Masters' candidate, Colorado School of Mines,
 *         Geophysics Department, 1994
 *
 * Implementation of "key=" option: Nils Maercklin,
 *         GeoForschungsZentrum (GFZ) Potsdam, Germany, 2002.
 *
 * References:
 *
 *    Embree, P.,1968, Diversity seismic record stacking method and system:
 *        U.S. patent 3,398,396.
 *    Gimlin, D. R., and Smith, J. W., 1980, A comparison of seismic trace 
 *        summing techniques: Geophysics, vol. 45, pages 1017-1041.
 *
 * Trace header fields accessed: ns, dt, key=keyword
 * Trace header fields modified: tracl
 */

/**************** end self doc ********************************/


segy intrace,outtrace;

void* main_sudivstack( void* args )
{   
    int nt;          /* number of samples input trace */
    char *key;       /* header key word from segy.h */
    char *type;      /* ... its type  */
    int index;       /* ... its index */
    Value val;       /* ... its value */
    int data=1;      /* flag no more data */
    int trnum;       /* current trace number */
    float winlen;    /* window length in seconds */
    int peak;        /* flag to use peak power */
    float entire;    /* flag to use entire trace as window */
    float dt;        /* sample rate in seconds */
    int ntwin;       /* number time samples in window */
    float totalpwa;  /* total power window a */
    float totalpwb;  /* total power window b */
    int j,i;         /* counters */
    float x;         /* variable */
    float avepwa;    /* average (or peak) power window a */
    float avepwb;    /* average (or peak) power window b */
    float maxval;    /* maximum data squared value in window */
    float *intpapw;  /* pointer array holding interpolated powers */
    int count;       /* counter */
    float *sumscale; /* pointer to sum of scalers array */
    float *sumdata;  /* pointer to summed data array */
    int nsegy;       /* number of bytes in segy */
    int newtracl=1;  /* tracl for stacked traces */
    float *xin;      /* pointer x interpolation */
    float *yin;      /* pointer y interpolation */
    float *xout;     /* pointer x interpolation */
    int z;           /* counter */
    
    /* Initialize */
    cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
    cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
    cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
    int argc = suArgs->argc;
    char **argv = suArgs->argv;
    cseis_su::csSUGetPars parObj;

    void* retPtr = NULL;  /*  Dummy pointer for return statement  */
    su2cs->setSUDoc( sdoc_sudivstack );
    if( su2cs->isDocRequestOnly() ) return retPtr;
    parObj.initargs(argc, argv);

    try {  /* Try-catch block encompassing the main function body */


    /* Get information from first trace */
    if(!(nsegy = cs2su->getTrace(&intrace)))
        throw cseis_geolib::csException("can not get first trace");

    nt = (int) intrace.ns;
    dt = ((double) intrace.dt)/1000000.0; /* microsecs to secs */ 
    entire = (float) ((nt * dt) - dt);    /* entire data length */
    
    /* Get parameters */
    if (!parObj.getparstring("key", &key)) key="tracf";
    if(!parObj.getparfloat("winlen", &winlen)) winlen = entire;
    if(!parObj.getparint("peak", &peak)) peak = 0;


    parObj.checkpars();

    ntwin = (int)(winlen/dt) + 1;
    if(ntwin > nt) ntwin = nt;            /* check data length */

    /* get key type and index */
    type = hdtype(key);
    index = getindex(key);

    /* allocate arrays */
    intpapw = ealloc1float(nt);
    sumscale = ealloc1float(nt);
    sumdata = ealloc1float(nt);
    xin=ealloc1float(2);
    yin=ealloc1float(2);
    xout=ealloc1float(ntwin);

    while(data) { 
        
        /* initialize arrays to zero */
        memset( (void *) sumscale, 0,nt * FSIZE);
        memset( (void *) sumdata, 0,nt * FSIZE);

        /* copy input header information to output header */
        memcpy((void*)&outtrace,(const void*)&intrace,nsegy);
    
        /* get current trace number (value of key from segy.h) */
        gethval(&intrace, index, &val);
        trnum = vtoi(type,val); 
    
        do {

            /* traces having same key=keyword are scaled and summed */
    
            /* calculate scaler in first window */
    
            if(peak) {              /* peak power */
                maxval = 0.0; 
                for(i=0;i<ntwin;++i) {
                    x = intrace.data[i];
                    maxval=((x*x)>maxval)?(x*x):maxval;
                }
                avepwa = maxval;
            } 
            else {                  /* average power */ 
                totalpwa = 0.0;      /*initialize */ 
                    for(i=0;i<ntwin;++i) {
                        x = intrace.data[i];
                        totalpwa += (x * x);
                    }
                    avepwa = totalpwa / ntwin;
                }
             
            /* interpolate from start to end of first window */
            for(i=0;i<ntwin;++i)
                intpapw[i] = avepwa;
    
            /* middle  windows */
            j = 1;
     
            while(ntwin * (j + 1) < nt) { 
                if(peak) {
                    maxval = 0.0;
                    for(i=ntwin * j;i<ntwin * (j+1);++i){
                        x = intrace.data[i];
                        maxval=((x*x)>maxval)?(x*x):maxval; 
                    }
                    avepwb = maxval;
                } 
                else { 
                    totalpwb = 0.0;
                    for(i=ntwin * j;i<ntwin * (j+1);++i) {
                        x = intrace.data[i];
                        totalpwb += (x * x);
                    }
                    avepwb = totalpwb / ntwin;
                }

                /* linear interpolate */ 
                /* between window a and window b */
                
                xin[0] = (ntwin * j) - 1;
                xin[1] = (ntwin * (j + 1)) - 1;
                yin[0] = avepwa;
                yin[1] = avepwb;
                z = 0;
                
                for(i=ntwin * j;i<ntwin * (j+1);++i) {
                    xout[z] = i; 
                    z++;
                }
                
                intlin(2,xin,yin,avepwa,avepwb,ntwin,xout,
                    &intpapw[ntwin * j]);

                /* update to next window */
                avepwa = avepwb;
                j++;
            }   /* end middle window while loop */
            
            /*  last window */
            if(ntwin < nt) {    
                if(peak) {
                    maxval = 0.0;
                    for(i=ntwin * j;i<nt;++i) {
                        x = intrace.data[i];
                        maxval = ((x * x) > maxval) ? 
                        (x * x) : maxval;
                    }
                    avepwb = maxval;
                } 
                else { 
                    count = 0;
                    totalpwb = 0.0;
                    for(i=ntwin * j;i<nt;++i) {
                        x = intrace.data[i];
                        totalpwb += (x * x);
                        count++;
                    }
                    avepwb = totalpwb / count;
                }
                /* interpolate last window */

                xin[0] = (ntwin * j) - 1;
                xin[1] = nt - 1;
                yin[0] = avepwa;
                yin[1] = avepwb;
                count = 0;
                z = 0;
                for(i=ntwin * j;i<nt;++i) {
                    xout[z] = i;
                    z++;
                    count++;
                }
                
                intlin(2,xin,yin,avepwa,avepwb,count,
                    xout,&intpapw[ntwin * j]);
            }    
            /* interpolation completed */
     
            /* multiply each data sample by its inverse ave power*/
            /* and sum if trace numbers the same                 */
            /* sum inverse ave powers for renormalization        */
            /* check not dividing by zero                        */
    
            for(i=0;i<nt;++i) {
                if(intpapw[i] != 0.0){
                    sumdata[i]+= intrace.data[i] / intpapw[i];
             
                    sumscale[i] += 1 / intpapw[i];
                }
            }    
            
            /* get next data trace and check trace number */
            if(!(nsegy = cs2su->getTrace(&intrace))) {
                data = 0;    /* set more data to false */
                break;       /* done no more traces  */
            }
            gethval(&intrace, index, &val);
        
        } while( vtoi(type,val) == trnum);
    
        /* renormalize data trace using inverse of sum of scalers */
        /* check not dividing by zero */
        for(i=0;i<nt;++i){
            if (sumscale[i] == 0.0) sumscale[i] = 1.0;

            outtrace.data[i] = sumdata[i] / sumscale[i];
        }

        outtrace.tracl = newtracl++; /* update header info */ 
        su2cs->putTrace(&outtrace);            /* write to disk */
    
    } /* end while continue loop */
    
    su2cs->setEOF();
    pthread_exit(NULL);
    return retPtr;

}
catch( cseis_geolib::csException& exc ) {
  su2cs->setError("%s",exc.getMessage());
  pthread_exit(NULL);
  return retPtr;
}
} /* end main loop */

} // END namespace
