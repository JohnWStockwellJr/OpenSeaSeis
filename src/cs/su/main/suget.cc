/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUGET: $Revision: 1.7 $ ; $Date: 2011/11/16 23:23:25 $                */

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
std::string sdoc_suget =
"									"
" SUGET  - Connect SU program to file descriptor for input stream.	"
"									"
"    suget fd=$1 | next_su_module					"
"									"
" This program is for interfacing \" outside processing systems \"	"
" with SU. Typically, an outside system would execute the su command file"
" and a file descriptor would be passed by an outside system to		"
" the su command file so that output data from the outside system	"
" could be piped into the su programs executing inside the command file."
"									"
" Example:    suget fd=$1 | next_su_module				"
"									"
"      fd=-1        file_descriptor_for_input_stream			"
"      verbose=0    minimal listing					"
"                   =1  asks for message with each trace processed.	"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suget {


/*
 * Author: John Anderson (visiting scholar from Mobil) July 1994
 */
/**************** end self doc ********************************/

int fgettrn(FILE *fp, segy *tp);

segy tr;
void* main_suget( void* args )
{
	FILE *fp;
	int fd,j=0,verbose;

	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suget );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	if(!parObj.getparint("fd",&fd)) fd=-1;
	if(!parObj.getparint("verbose",&verbose)) verbose=0; 
        parObj.checkpars();
	warn("File descriptor passed to suget = %d",fd);
	if( (fp = (FILE *) fdopen(fd,"r"))==NULL) throw cseis_geolib::csException("Bad file descriptor");
	warn("About to read first trace");
	if(!fgettr(fp,&tr)) throw cseis_geolib::csException("Can't get first trace");

	do{
		if(verbose>0){
			warn("read trace %d",j);
			j++;
		}

		su2cs->putTrace(&tr);

	} while(fgettr(fp,&tr));

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
