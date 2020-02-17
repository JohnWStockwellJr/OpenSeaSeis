/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUMEAN: $Revision: 1.5 $ ; $Date: 2011/11/16 17:24:58 $	*/

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

/*************************** self documentation ******************************/
std::string sdoc_sumean =
"									"
" SUMEAN - get the mean values of data traces				"	
"									"
" sumean < stdin > stdout [optional parameters] 			"
"									"
" Required parameters:							"
"   power = 2.0		mean to the power				"
"			(e.g. = 1.0 mean amplitude, = 2.0 mean energy)	"
"									"
" Optional parameters: 							"
"   verbose = 0		writes mean value of section to outpar	   	"
"			= 1 writes mean value of each trace / section to"
"				outpar					"
"   outpar=/dev/tty   output parameter file				"
"   abs = 1             average absolute value "
"                       = 0 preserve sign if power=1.0"
"			 						"
" Notes:			 					"
" Each sample is raised to the requested power, and the sum of all those"
" values is averaged for each trace (verbose=1) and the section.	"
" The values power=1.0 and power=2.0 are physical, however other powers	"
" represent other mathematical L-p norms and may be of use, as well.	"
"			 						"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sumean {


/* Credits:
 *  Bjoern E. Rommel, IKU, Petroleumsforskning / October 1997
 *		    bjorn.rommel@iku.sintef.no
 */

/**************** end self doc ***********************************************/


/* Globals */
segy tr;

void* main_sumean( void* args )
{
  char *outpar ="";	   /* name of output file			*/
  FILE *outparfp = NULL;   /* file pointer of output file		*/
  int verbose = 0;	   /* flag for printing extra information	*/
  int itr = 0;   	   /* trace number				*/
  int it = 0;		   /* sample number				*/
  float power = 0.0;       /* mean to the power of			*/
  float tmean = 0.0;	   /* average mean of trace			*/
  float gmean = 0.0;	   /* average mean of section			*/
  int abs = 1;             /* absolute value flag */


  /* Initialize */
  cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
  cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
  cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
  int argc = suArgs->argc;
  char **argv = suArgs->argv;
  cseis_su::csSUGetPars parObj;

  void* retPtr = NULL;  /*  Dummy pointer for return statement  */
  su2cs->setSUDoc( sdoc_sumean );
  if( su2cs->isDocRequestOnly() ) return retPtr;
  parObj.initargs(argc, argv);

  try {  /* Try-catch block encompassing the main function body */


  /* Get optional parameters */
  if (!parObj.getparint ("verbose", &verbose))   verbose = 0;
  if ((verbose < 0) || (verbose > 1))   throw cseis_geolib::csException("unknown mode for verbose!");
  if (!parObj.getparstring ("outpar", &outpar))   outpar = "/dev/tty" ;
  outparfp = efopen (outpar, "w");
  if (!parObj.getparint ("abs", &abs))   abs = 0;
  if (!parObj.getparfloat ("power", &power))   power = 2.0;

  parObj.checkpars();
  /* Get info from first trace */
  if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");

  /* Initialize mean value of section */
  gmean = 0.0;

  /* Loop through traces */
  itr = 0;
  do {
    /* Reset */
    itr++ ;	/* count traces */
    tmean = 0.0;   /* reset mean value of trace */

    /* Loop through samples */
    for (it = 0; it < tr.ns; ++it) {
    
	/* Raise sample to the requested power, add to mean value of trace */
        if( abs ){

	   tmean += pow (fabs (tr.data[it]), power);

        }else{

	   tmean += pow ((tr.data[it]), power);

        }
    
    }
    
    /* Average mean value of trace */
    tmean = pow (tmean / tr.ns, 1.0 / power);

    /* Add to the mean value of section */
    gmean += tmean;

    /* Print mean value of trace */
    if (verbose == 1)   
	fprintf (outparfp, "trace: %i   mean value: %e\n", itr, tmean);

  } while (cs2su->getTrace(&tr)); 

  /* Average mean value of section */
  gmean = gmean / itr;

  /* Print mean value of section */
  fprintf (outparfp, "global mean: %e\n", gmean);
  
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
