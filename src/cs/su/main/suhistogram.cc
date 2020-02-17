/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUHISTOGRAM: $Revision: 1.10 $ ; $Date: 2011/11/16 17:24:58 $	*/

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
#include "math.h" 

/*********************** self documentation **********************/
std::string sdoc_suhistogram =
"								"
" SUHISTOGRAM - create histogram of input amplitudes		"
"								"
"    suhistogram <in.su >out.dat				"
"								"
" Required parameters:						"
" min=		minimum bin 					"
" max=		maximum bin 					"
" bins=		number of bins					"
"								"
" Optional parameters						"
" trend=0	=0 1-D histogram				"
"	   =1 2-D histogram picks on cumulate			"
"	   =2 2-D histogram in trace format			"
"								"
" clip=     threshold value to drop outliers			"
"								"
" dt=	sample rate in feet or milliseconds.  Defaults  to	"
"    	tr.dt*1e-3					  	"
" datum=  header key to get datum shift if desired (e.g. to	"
"	 hang from water bottom)			    	"
"								"
" Notes:							"
" trend=0 produces a two column ASCII output for use w/ gnuplot."
" Extreme values are counted in the end bins.			"
"								"
" trend=1 produces a 6 column ASCII output for use w/ gnuplot   "
" The columns are time/depth and picks on the cumulate		"
" at 2.28%, 15.87%, 50%, 84.13% & 97.72% of the total points    "
" corresponding to the median and +- 1 or 2 standard deviations "
" for a Gaussian distribution.					"
"								"
" trend=2 produces an SU trace panel w/ one trace per bin that  "
" can be displayed w/ suximage, etc.				"
"							   	"
" Example for plotting with xgraph:				"
" suhistogram < data.su min=MIN max=MAX bins=BINS |		"
" a2b n1=2 | xgraph n=BINS nplot=1			 	"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suhistogram {


/*
 * Author: Reginald H. Beardsley  2006   rhb@acm.org
 * 
 */

/**************** end self doc ********************************/

segy     tr;

void* main_suhistogram( void* args )
{

   int i;  
   int j;  
   int k;  
   int bins;
   int trend=0;

   int p1;
   int p2;
   int p3;
   int p4;
   int p5;

   float m;
   int sum;
   float *a=NULL;
   float *b=NULL;
   float c;

   float *cnt=NULL;
   int **histogram=NULL;

   float min;
   float max;
   float bin;
   Value val;
   float fval;
   float dt;

   float clip=1e38;

   cwp_String datum=0;

   /* Initialize */
   cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
   cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
   cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
   int argc = suArgs->argc;
   char **argv = suArgs->argv;
   cseis_su::csSUGetPars parObj;

   void* retPtr = NULL;  /*  Dummy pointer for return statement  */
   su2cs->setSUDoc( sdoc_suhistogram );
   if( su2cs->isDocRequestOnly() ) return retPtr;
   parObj.initargs(argc, argv);

   try {  /* Try-catch block encompassing the main function body */


   /* Get info from first trace */
   if( !cs2su->getTrace(&tr) ){
      throw cseis_geolib::csException("Can't get first trace \n");
   }

   /* Get parameters */
   if( !parObj.getparfloat("min", &min) ){
      throw cseis_geolib::csException("min must be specified\n");
   }

   if( !parObj.getparfloat("max", &max) ){
      throw cseis_geolib::csException("max must be specified\n");
   }

   if( !parObj.getparint("bins", &bins) ){
      throw cseis_geolib::csException("bins must be specified\n");
   }

   parObj.getparint("trend", &trend);

   parObj.getparstring( "datum" ,&datum );

   parObj.getparfloat( "clip" ,&clip );

   if( !parObj.getparfloat( "dt" ,&dt ) ){
      dt = tr.dt*1e-3;
   }

   parObj.checkpars();
   bin = (max-min) / bins;

   cnt = ealloc1float( tr.ns );
   a = ealloc1float( tr.ns );
   b = ealloc1float( tr.ns );
   memset( cnt ,0 ,tr.ns*sizeof(float) );
   memset( a   ,0 ,tr.ns*sizeof(float) );
   memset( b   ,0 ,tr.ns*sizeof(float) );

   if( trend ){

      /*-------------------*/
      /* allocate a buffer */
      /*-------------------*/

      histogram=ealloc2int( bins ,tr.ns );

      for( i=0; i<tr.ns; i++ ){
	 memset( histogram[i] ,0 ,bins*sizeof(int) );
      }

      /*----------------------*/
      /* accumulate histogram */
      /*----------------------*/

      do{

	 if( datum ){
	    gethdval(&tr ,datum ,&val);
	    fval = vtof(hdtype(datum), val);

	 }else{
	    fval = 0.0;

	 }

	 for( i=0; i<tr.ns; i++ ){
   
	    if( tr.data[i] < clip ){

	       j = (tr.data[i] - min) / bin;
	       k = (i*dt - fval)/dt;
      
	       if( j >= 0 && j < bins 
		&& k >= 0 && k < tr.ns ){
		  histogram[k][j]++;
		  a[k] += tr.data[i]*tr.data[i];
		  b[k] += tr.data[i];
		  cnt[k]++;
	       } 
	    }
	 }
   
      }while( cs2su->getTrace(&tr) );

      if( trend == 1 ){

	 /*-----------------*/
	 /* pick key points */
	 /*-----------------*/

	 for( i=0; i<tr.ns; i++ ){


	    if( cnt[i] > 0 ){

	       p1 = 0;
	       p2 = 0;
	       p3 = 0;
	       p4 = 0;
	       p5 = 0;

	       sum = 0;

	       printf( "%f " ,i*tr.dt/1000.0 );

	       for( j=0; j<bins; j++ ){

		  sum += histogram[i][j];

		  m = (100.0*sum ) / cnt[i]; 
		     
		  if( !p1 && m >= 2.28 ){
		     p1 = 1;
		     printf( "%f " ,min + j*bin );

		  }else if( !p2 && m >= 15.87 ){
		     p2 = 1;
		     printf( "%f " ,min + j*bin );

		  }else if( !p3 && m >= 50.00 ){
		     p3 = 1;
		     printf( "%f " ,min + j*bin );

		  }else if( !p4 && m >= 84.13 ){
		     p4 = 1;
		     printf( "%f " ,min + j*bin );

		  }else if( !p5 && m >= 97.72 ){
		     p5 = 1;
		     printf( "%f " ,min + j*bin );
		  }

	       }

	       c = sqrt( (cnt[i]*a[i] - b[i]*b[i])/(cnt[i]*(cnt[i]-1)));

	       printf( "%f " ,c );

	       printf( "\n" );

	    }

	 }

      }else{

	 for( j=0; j<bins; j++ ){

	    for( i=0; i<tr.ns; i++ ){

	       if( cnt[i] > 0 ){
		  tr.data[i] = histogram[i][j] / (float) cnt[i];

	       }else{
		  tr.data[i] = histogram[i][j];

	       }
	       
	    }
	    tr.tracl = j;
	    cs2su->putTrace(&tr );

	  }


      }

   }else{

      histogram=ealloc2int( 1 ,bins );
      memset( histogram[0] ,0 ,bins*sizeof(int) );
   
      /* Loop over traces */
      do{
	 for( i=0; i<tr.ns; i++ ){
   
	    j = (tr.data[i] - min) / bin;
	    j = j < 0 ? 0 : j;
	    j = j > bins-1 ? bins-1 : j;
   
	    histogram[0][j]++;
	    cnt[0]++;
   
	 }
   
      }while( cs2su->getTrace(&tr) );

      for( i=0; i<bins; i++ ){
   
	 printf( "%15f " ,min+i*bin		);
	 printf( "%15f " , (float)  histogram[0][i]/cnt[0] );
	 printf( "\n" );
      } 



   }

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
