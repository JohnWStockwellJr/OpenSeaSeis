/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* : $Revision: 1.5 $ ; $Date: 2011/11/16 17:43:20 $	*/

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
#include "par.h"
#include "su.h"

/*********************** self documentation ***************************/

std::string sdoc_segyhdrmod =
" SEGYHDRMOD - replace the text header on a SEGY file		"
"								"
"   segyhdrmod text=file data=file				"
"								"
"   Required parameters:					"
"								"
"   text=      name of file containing new 3200 byte text header"
"   data=      name of file containing SEGY data set		"
"								"
" Notes:							"
" This program simply does a replacement of the content of the first"
" 3200 bytes of the SEGY file with the contents of the file specified"
" by the text= parameter. If the text header in the SEGY standard"
" ebcdic format, the user will need to supply an ebcdic format file"
" as the text=  as input file. A text file may be converted from"
" ascii to ebcdic via:						"
"   dd if=ascii_filename of=ebcdic_filename conv=ebcdic ibs=3200 count=1"
" or from ebcdic to ascii via:					"
"   dd if=ebcdic_filename of=ascii_filename ibs=3200 conv=ascii count=1"
" 								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace segyhdrmod {



/*====================================================================*\

   sgyhdrmod - replace the text header on a SEGY data file in place

   This program only reads and writes 3200 bytes

   Reginald H. Beardsley                            rhb@acm.org

\*====================================================================*/
/************************** end self doc ******************************/
      

void* main_segyhdrmod( void* args )
{

   FILE *txtfp=NULL;   /* file pointer for new text */
   FILE *datfp=NULL;   /* file pointer for data file */

   int n;

   char *text=NULL;
   char *data=NULL;

   char buf[3200];

   /*------------*/
   /* Initialize */
   /*------------*/

   cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
   cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
   cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
   int argc = suArgs->argc;
   char **argv = suArgs->argv;
   cseis_su::csSUGetPars parObj;

   void* retPtr = NULL;  /*  Dummy pointer for return statement  */
   su2cs->setSUDoc( sdoc_segyhdrmod );
   if( su2cs->isDocRequestOnly() ) return retPtr;
   parObj.initargs(argc, argv);

   try {  /* Try-catch block encompassing the main function body */


   if( !parObj.getparstring( "text" ,&text ) ){
      throw cseis_geolib::csException( "missing text header filename" );
   }
 
   if( !parObj.getparstring( "data" ,&data ) ){
      throw cseis_geolib::csException( "missing data filename" );
   }
   parObj.checkpars();
 
   /*------------------*/
   /* Open input files */
   /*------------------*/

   if( !(txtfp = fopen(text, "rb")) ){
      throw cseis_geolib::csException( "unable to open %s" ,text );
   }

   if( !(datfp = fopen(data, "rb+")) ){
      throw cseis_geolib::csException( "unable to open %s" ,data );
   }

   /*---------------------------*/
   /* rewrite text header block */
   /*---------------------------*/

   if( (n=fread( buf ,1 ,sizeof(buf) ,txtfp )) != sizeof(buf) ){
      throw cseis_geolib::csException( "unable to read new text header" );
   }

   if( (n=fwrite( buf ,1 ,sizeof(buf) ,datfp )) != sizeof(buf) ){
      throw cseis_geolib::csException( "write of new text header failed!!!!!" );
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
