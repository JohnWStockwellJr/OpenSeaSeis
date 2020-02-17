/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUCHART: $Revision: 1.19 $ ; $Date: 2011/11/16 22:10:29 $	*/

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
std::string sdoc_suchart =
" 								"
" SUCHART - prepare data for x vs. y plot			"
" 								"
" suchart <stdin >stdout key1=sx key2=gx			"
" 								"
" Required parameters:						"
" 	none							"
" 								"
" Optional parameters:						"
" 	key1=sx  	abscissa 				"
" 	key2=gx		ordinate				"
"	outpar=null	name of parameter file			"
" 								"
" The output is the (x, y) pairs of binary floats		"
" 								"
" Examples:							"
" suchart < sudata outpar=pfile >plot_data			"
" psgraph <plot_data par=pfile title=\"CMG\" \\			"
"	linewidth=0 marksize=2 mark=8 | ...			"
" rm plot_data 							"
" 								"
" suchart < sudata | psgraph n=1024 d1=.004 \\			"
"	linewidth=0 marksize=2 mark=8 | ...			"
" 								"
" fold chart: 							"
" suchart < stacked_data key1=cdp key2=nhs |			"
"            psgraph n=NUMBER_OF_TRACES d1=.004 \\		"
"	linewidth=0 marksize=2 mark=8 > chart.ps		"
" 								"
" 								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suchart {


/* Credits:
 *	SEP: Einar Kjartansson
 *	CWP: Jack K. Cohen
 *
 * Notes:
 *	The vtof routine from valpkge converts values to floats.
 */
/**************** end self doc ***********************************/


segy tr;

void* main_suchart( void* args )
{
	cwp_String key1,  key2;	/* x and y key header words	*/
	Value  val1,  val2;	/* ... their values		*/
	cwp_String type1, type2;/* ... their types		*/
	int index1, index2;	/* ... their indices in hdr.h	*/
	float x, y;		/* temps to hold current x & y 	*/
	cwp_String outpar;	/* name of par file		*/
	register int npairs;	/* number of pairs found	*/


	/* Hook up getpars */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suchart );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Prevent byte codes from spilling to screen */
	if (isatty(STDOUT)) throw cseis_geolib::csException("must redirect or pipe binary output");


	/* Get parameters */
	if (!parObj.getparstring("key1", &key1))	key1 = "sx";
	if (!parObj.getparstring("key2", &key2))	key2 = "gx";

	type1 = hdtype(key1);
	type2 = hdtype(key2);

	index1 = getindex(key1);
	index2 = getindex(key2);


	/* Loop over traces */
	npairs = 0;
	while(cs2su->getTrace(&tr)) {

		gethval(&tr, index1, &val1);
		gethval(&tr, index2, &val2);

		x = vtof(type1, val1);
		y = vtof(type2, val2);

		efwrite(&x, FSIZE, 1, stdout);
		efwrite(&y, FSIZE, 1, stdout);

		++npairs;
	}


	/* Make parfile if needed */
	if (parObj.getparstring("outpar", &outpar))
		fprintf(efopen(outpar, "w"),
			"n=%d label1=%s label2=%s\n",
			npairs, key1, key2);

        parObj.checkpars();

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
