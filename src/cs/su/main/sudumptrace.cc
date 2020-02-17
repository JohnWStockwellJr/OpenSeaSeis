/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUDUMPTRACE: $Revision: 1.3 $ ; $Date: 2011/11/16 22:10:29 $	*/

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
  #include <signal.h>

}

/*********************** self documentation **********************/
std::string sdoc_sudumptrace =
"                                                                   "
" SUDUMPTRACE - print selected header values and data.              "
"               Print first num traces.                             "
"               Use SUWIND to skip traces.                          "
"                                                                   "
" sudumptrace < stdin [> ascii_file]                                "
"                                                                   "
" Optional parameters:                                              "
"     num=4                    number of traces to dump             "
"     key=key1,key2,...        key(s) to print above trace values   "
"     hpf=0                    header print format is float         "
"                              =1 print format is exponential       "
"                                                                   "
" Examples:                                                         "
"   sudumptrace < inseis.su            PRINTS: 4 traces, no headers "
"   sudumptrace < inseis.su key=tracf,offset                        "
"   sudumptrace < inseis.su num=7 key=tracf,offset > info.txt       "
"   sudumptrace < inseis.su num=7 key=tracf,offset hpf=1 > info.txt "
"                                                                   "
" Related programs: suascii, sugethw                                "
"                                                                   "
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sudumptrace {


/* Credits:
 *   MTU: David Forel, Jan 2005
 *
 * Trace header field accessed: nt, dt, delrt
 */
/**************** end self doc ***********************************/

/* subroutine prototypes */
void dump(float *data, float dt, float *hedr, cwp_String *key,
          float delrt, int nkeys, int ntr, int nt, int hpf) ;
static void closefiles(void) ;

/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ] ;   /* filename for the file of traces */
char headerfile[BUFSIZ] ;  /* filename for the file of headers */
FILE *tracefp ;            /* fp for trace storage file */
FILE *headerfp ;           /* fp for header storage file */

segy tr ;

void* main_sudumptrace( void* args )
{
   char *tmpdir ;                 /* directory path for tmp files */
   cwp_Bool istmpdir=cwp_false ;  /* true for user given path */
   float *hedr ;                  /* the headers */
   float *data ;                  /* the data */

   int nt ;                       /* number of trace samples */
   float dt ;                     /* sample interval, sec */
   float delrt ;                  /* delay recording time, sec */
   cwp_String key[SU_NKEYS] ;     /* array of keywords */
   cwp_String type ;              /* key string type */
   int nkeys ;                    /* number of keywords */
   int ikey,ntr = 0 ;	          /* counters */
   int num ;                      /* number of traces to dump */
   int numtr = 4 ;                /* number of traces to dump */
   int hpf ;                      /* header print format */

   /* Initialize */
   cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
   cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
   cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
   int argc = suArgs->argc;
   char **argv = suArgs->argv;
   cseis_su::csSUGetPars parObj;

   void* retPtr = NULL;  /*  Dummy pointer for return statement  */
   su2cs->setSUDoc( sdoc_sudumptrace );
   if( su2cs->isDocRequestOnly() ) return retPtr;
   parObj.initargs(argc, argv);

   try {  /* Try-catch block encompassing the main function body */


   /* Look for user-supplied tmpdir */
   if (!parObj.getparstring("tmpdir",&tmpdir) &&
       !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
   if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
       throw cseis_geolib::csException("you can't write in %s (or it doesn't exist)", tmpdir);

   /* Get values from first trace */
   if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
   nt = (int) tr.ns ;                  /* Get nt */
   dt = ((double) tr.dt)/1000000.0 ;   /* microsecs to secs */
   if (!dt) parObj.getparfloat("dt", &dt) ;
   if (!dt) CSMUSTGETPARFLOAT("dt", &dt) ;
   delrt = ((double) tr.delrt)/1000.0 ; /* millisecs to secs */

   /* Get parameters */
   if (parObj.getparint ("num", &num)) numtr = num ;
   if ((nkeys=parObj.countparval("key"))!=0) parObj.getparstringarray("key",key) ;
   hedr = ealloc1float(nkeys*numtr) ;  /* make space for headers */
   if (!parObj.getparint ("hpf", &hpf)) hpf = 0 ;
   parObj.checkpars();

   /* Store traces, headers in tempfiles */
   if (STREQ(tmpdir,""))
   {
      tracefp = etmpfile();
      headerfp = etmpfile();

      do
      {
         ++ntr;
         efwrite(&tr, HDRBYTES, 1, headerfp);
         efwrite(tr.data, FSIZE, nt, tracefp);

         /* Get header values */
         for (ikey=0; ikey<nkeys; ++ikey)
         {
            Value val;
            float fval;

            gethdval(&tr, key[ikey], &val) ;
            type = hdtype(key[ikey]) ;
            fval = vtof(type,val) ;
            hedr[(ntr-1)*nkeys+ikey] = fval ;
         }

      }
      while (ntr<numtr  &&  cs2su->getTrace(&tr)) ;

   }
   else  /* user-supplied tmpdir */
   {
      char directory[BUFSIZ];
      strcpy(directory, tmpdir);
      strcpy(tracefile, temporary_filename(directory));
      strcpy(headerfile, temporary_filename(directory));
      /* Handle user interrupts */
      signal(SIGINT, (void (*) (int)) closefiles);
      signal(SIGQUIT, (void (*) (int)) closefiles);
      signal(SIGHUP,  (void (*) (int)) closefiles);
      signal(SIGTERM, (void (*) (int)) closefiles);
      tracefp = efopen(tracefile, "w+");
      headerfp = efopen(headerfile, "w+");
      istmpdir=cwp_true;      

      do
      {
         ++ntr;
         efwrite(&tr, HDRBYTES, 1, headerfp);
         efwrite(tr.data, FSIZE, nt, tracefp);

         /* Get header values */
         for (ikey=0; ikey<nkeys; ++ikey)
         {
            Value val;
            float fval;

            gethdval(&tr, key[ikey], &val) ;
            type = hdtype(key[ikey]) ;
            fval = vtof(type,val) ;
            hedr[(ntr-1)*nkeys+ikey] = fval ;
         }

      }
      while (ntr<numtr  &&  cs2su->getTrace(&tr)) ;

   }

   /* Rewind after read, allocate space */
   erewind(tracefp);
   erewind(headerfp);
   data = ealloc1float(nt*ntr);

   /* Load traces into data and close tmpfile */
   efread(data, FSIZE, nt*ntr, tracefp);
   efclose(tracefp);
   if (istmpdir) eremove(tracefile);

   rewind(headerfp);
   rewind(tracefp);

   /* Do trace work */
   dump(data, dt, hedr, key, delrt, nkeys, ntr, nt, hpf) ;

   /* close */
   efclose(headerfp);
   if (istmpdir) eremove(headerfile);

   free1(hedr) ;
   free1(data) ;

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


void dump(float *data, float dt, float *hedr, cwp_String *key,
          float delrt, int nkeys, int ntr, int ntime, int hpf)
/*
  Dump headers and traces in column format; one trace per column.
  INPUT:
    data   array of trace data values
    dt     trace sample interval
    delrt  delay recording time, sec
    hedr   array of trace headers
    key    array of key names
    nkeys  number of keys (headers) to print above trace values
    ntr    number of traces to dump
    ntime  number of time samples on trace
    hpf    header print format flag: 0=float, 1=exponential
  OUTPUT:  none
*/
{
   int i,j,k, m ;         /* counters */

   printf("\nnum traces = %d    num samples = %d \n\n",ntr,ntime) ;

   for (k=0; k<nkeys; ++k)                     /* Print headers */
   {
      printf("%17s  ", key[k]) ;
      for (m=0; m<ntr; ++m)
      {
         if (hpf==0)
         {
            printf("%11.4f\t", hedr[m*nkeys + k]) ;
         }
         else
         {
            printf("%11.4e\t", hedr[m*nkeys + k]) ;
         }
      }
      putchar('\n') ;
   }

   putchar('\n') ;

   printf("\nCounter    Time     Values\n") ;  /* Column titles */

   for (i=1; i<=ntime; ++i)                    /* Print trace values */
   {
      printf(" %6d ", i) ;
      printf(" %8.3f   ", dt*(i)+delrt) ;
      for (j=1; j<=ntr; ++j)
      {
         printf("%11.4e\t", data[(j-1)*ntime+(i-1)]) ;
      }
      putchar('\n') ;
   }
   putchar('\n') ;

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
