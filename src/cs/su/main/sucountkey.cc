/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUCOUNTKEY: $Revision: 1.6 $ ; $Date: 2011/11/16 22:10:29 $	*/

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


/*********************** self documentation ******************************/
std::string sdoc_sucountkey =
"									"
" SUCOUNTKEY - COUNT the number of unique values for a given KEYword.	"
"									"
" sucountkey < input.su key=[sx,gx,cdp,...]				"
" Required parameter:							"
" key=			array of SU header keywords being counted	"
" Optional parameters:							"
" verbose=1		chatty, =0 just print keyword number		"
" Example:								"
"	  suplane | sucountkey key=tracl,tracr,offset			"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sucountkey {


/*
 * Credits: Baoniu Han, bhan@mines.edu, Nov, 2000
 */
/**************** end self doc ********************************/

/* Globals variables*/
segy tr;

/* internal structure */
typedef struct {
	int num;
	Value *val;
} my_Count;

void* main_sucountkey( void* args )
{
	cwp_String key[SU_NKEYS];	/* array of keywords		*/
	cwp_String type[SU_NKEYS];	/* array of keywords		*/
	int index[SU_NKEYS];		/* name of type of getparred key*/ 
	my_Count count_table[SU_NKEYS];	/* count table for keywords	*/

	int verbose;			/* verbose flag			*/

	int ikey;			/* key counter			*/
	int nkeys;			/* number of header fields set	*/
	int i;				/* count index			*/
	Value val0;			/* value of key field		*/
	long buf_size;
	
	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sucountkey );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

	
	/* Get "key" values */ 
	if ((nkeys=parObj.countparval("key"))!=0) {
		parObj.getparstringarray("key",key); 
	} else {
		key[0]="cdp"; 
	}
	if (!parObj.getparint("verbose", &verbose))	verbose = 1;

        parObj.checkpars();
	/* you can increase this number to change the buffer size */
	/* for each key word */
	buf_size=100000L;
	
	/* get types and indexes corresponding to the keys */
	for (ikey=0; ikey<nkeys; ++ikey) {
		type[ikey]=hdtype(key[ikey]); 
		index[ikey]=getindex(key[ikey]);
		count_table[ikey].num=0;
		count_table[ikey].val=(Value *) malloc(sizeof(Value)*buf_size);
	}


	/* get info from first trace */
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace"); 
	for(ikey=0;ikey<nkeys;++ikey){
		count_table[ikey].num=1;
		gethval(&tr, index[ikey], &count_table[ikey].val[0]);
	}



	/* loop over traces */
	while (cs2su->getTrace(&tr)) {
		
	  for (ikey=0; ikey<nkeys; ++ikey) {
		gethval(&tr,index[ikey],&val0); 
		for(i=0;i<count_table[ikey].num;i++)
	 	    if (valcmp(type[ikey], val0, count_table[ikey].val[i])==0) break;

	  if (i==count_table[ikey].num){				
		gethval(&tr, index[ikey], &count_table[ikey].val[count_table[ikey].num]);		
		count_table[ikey].num++;
	  } 

	}

	}

	if (verbose==1)
		warn("The counting of key word is finished, the total number for specified header is as follow:");

	for (ikey=0;ikey<nkeys;++ikey){
		if (verbose==1)
			warn(" key %s has %d unique entries in this data", key[ikey], count_table[ikey].num);
		else
			fprintf(stdout,"%s %d ", key[ikey], count_table[ikey].num);

		
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
