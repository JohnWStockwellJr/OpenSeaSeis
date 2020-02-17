/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUKEYCOUNT: $Revision: 1.3 $ ; $Date: 2011/11/16 22:10:29 $	*/

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
std::string sdoc_sukeycount =
"                                                             "
" SUKEYCOUNT - sukeycount writes a count of a selected key    "
"                                                             "
"   sukeycount key=keyword < infile [> outfile]                  "
"                                                             "
" Required parameters:                                        "
" key=keyword      One key word.                                 "
"                                                             "
" Optional parameters:                                        "
" verbose=0  quiet                                            "
"        =1  chatty                                           "
"                                                             "
" Writes the key and the count to the terminal or a text      "
"   file when a change of key occurs. This does not provide   "
"   a unique key count (see SUCOUNTKEY for that).             "
" Note that for key values  1 2 3 4 2 5                       "
"   value 2 is counted once per occurrence since this program "
"   only recognizes a change of key, not total occurrence.    "
"                                                             "
" Examples:                                                   "
"    sukeycount < stdin key=fldr                              "
"    sukeycount < stdin key=fldr > out.txt                    "
"                                                             "
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sukeycount {


/* Credits:
 *
 *   MTU: David Forel, Jan 2005
 */
/**************** end self doc ***********************************/

segy tr ;

void* main_sukeycount( void* args )
{
   cwp_String key[SU_NKEYS] ;  /* array of keywords */
   int nkeys ;            /* number of keywords to retrieve */
   int iarg ;             /* arguments in argv loop */
   int countkey = 0 ;     /* counter of keywords in argc loop */
   int verbose = 0 ;      /* verbose ? */
   cwp_String type1=NULL;     /* key string */
   float sort, sortold ;  /* for comparing new/old key values */
   int isort, isortold ;  /* for comparing new/old key values */
   int gatherkount ;      /* counter of traces within key */
   int itotal = 0 ;       /* counter of total traces */

   /* Initialize */
   cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
   cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
   cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
   int argc = suArgs->argc;
   char **argv = suArgs->argv;
   cseis_su::csSUGetPars parObj;

   void* retPtr = NULL;  /*  Dummy pointer for return statement  */
   su2cs->setSUDoc( sdoc_sukeycount );
   if( su2cs->isDocRequestOnly() ) return retPtr;
   parObj.initargs(argc, argv);

   try {  /* Try-catch block encompassing the main function body */


   sortold  = -99999. ;
   isortold  = -99999 ;

   /* Get key values */
   if (!parObj.getparint("verbose",&verbose)) verbose=0 ;
   if ((nkeys=parObj.countparval("key"))!=0) {
      parObj.getparstringarray("key",key) ;
   } else {
      /* support old fashioned method for inputting key fields */
      /* as single arguments:  sukeycount key1 */
      if (argc==1) throw cseis_geolib::csException("must set one key value!") ;

      for (iarg = 1; iarg < argc; ++iarg)
      {
         cwp_String keyword ;  /* keyword */

         keyword = argv[iarg] ;

         if (verbose) warn("argv=%s",argv[iarg]);
         /* get array of types and indexes to be set */
         if ((strncmp(keyword,"output=",7)!=0)) {
            key[countkey] = keyword ;
            ++countkey ;
         }
         if (countkey==0) throw cseis_geolib::csException("must set one key value!") ;
         if (countkey>1) throw cseis_geolib::csException("must set only one key value!") ;
      }
      nkeys=countkey;
   }
   parObj.checkpars();
   if (nkeys>1) throw cseis_geolib::csException("must set only one key value!") ;

   printf("\n") ;

   /* Loop over traces */
   gatherkount = 0 ;
   while (cs2su->getTrace(&tr)) {
      /* Do not loop over keys because only one is used */
      Value vsort ;
      gethdval(&tr, key[0], &vsort) ;
      type1 = hdtype(key[0]) ;

      if (*type1 == 'f')  /* float header */ {
         sort  = vtof(type1,vsort) ;

         /* Don't write just because first trace is new */
         if ( itotal == 0 ) sortold = sort ;

         if ( sort != sortold ) {
            printf(" %8s = %f", key[0], sortold) ;
            printf("      has  %d  trace(s)", gatherkount) ;
            printf("\n") ;
            sortold = sort ;
            gatherkount = 0 ;
         }
         ++gatherkount ;
         ++itotal ;
      } else  /* non-float header */ {
         isort = vtoi(type1,vsort) ;

         /* Don't write just because first trace is new */
         if ( itotal == 0 ) isortold = isort ;

         if ( isort != isortold )
         {
            printf(" %8s = %d", key[0], isortold) ;
            printf("      has  %d  trace(s)", gatherkount) ;
            printf("\n") ;
            isortold = isort ;
            gatherkount = 0 ;
         }
         ++gatherkount ;
         ++itotal ;
      }
   }

   /* Write after last trace is read */
   if (*type1 == 'f')
   {
      printf(" %8s = %f", key[0], sortold) ;
      printf("      has  %d  trace(s)", gatherkount) ;
      printf("\n") ;
   }
   else
   {
      printf(" %8s = %d", key[0], isortold) ;
      printf("      has  %d  trace(s)", gatherkount) ;
      printf("\n") ;
   }
   printf("\n        %d trace(s) read\n\n", itotal) ;

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
