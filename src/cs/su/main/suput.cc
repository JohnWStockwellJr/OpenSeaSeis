/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUPUT: $Revision: 1.7 $ ; $Date: 2011/11/16 23:23:25 $                */

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
std::string sdoc_suput =
"									"
" SUPUT - Connect SU program to file descriptor for output stream.	"
"									"
"       su_module | suput fp=$1						"
"									"
" This program is for interfacing \" outside processing systems \"	"
" with SU. Typically, the outside system would execute the SU command file."
" The outside system provides the file descriptor it would like to read	"
" from to the command file to be an argument for suput.			"
"									"
" Example: su_module | suput fp=$1					"
"									"
"       fd=-1       file_descriptor_for_output_stream_from_su		"
"       verbose=0   minimal listing					"
"                   =1  asks for message with each trace processed.	"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suput {


/*
 * Author: John Anderson (visiting scholar from Mobil) July 1994
 */

/**************** end self doc ********************************/

void fputtrn(FILE *fp, segy *tp);

segy tr;
void* main_suput( void* args )
{
	FILE *fp;
	int fd,nread=0,verbose;

	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suput );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	if(!parObj.getparint("fd",&fd)) fd=-1;
	fprintf(stderr,"File descriptor passed to suput = %d\n",fd);

	if(!parObj.getparint("verbose",&verbose)) verbose=0; 
        parObj.checkpars();

	if( (fp = (FILE *) fdopen(fd,"w"))==NULL) throw cseis_geolib::csException("Bad file descriptor \n");

	if(!cs2su->getTrace(&tr)) throw cseis_geolib::csException("Can't get first trace \n");
	do{
		if(verbose>0){
			warn("suput: read input traces %d", nread);
			nread++;
		}

		fputtr(fp,&tr);
	} while(cs2su->getTrace(&tr));

	(void) fclose(fp);
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
