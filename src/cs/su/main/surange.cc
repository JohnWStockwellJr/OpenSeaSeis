/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SURANGE: $Revision: 1.18 $ ; $Date: 2011/11/16 22:10:29 $  */

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
std::string sdoc_surange =
" 								"
" SURANGE - get max and min values for non-zero header entries	"
" 								"
" surange <stdin	 					"
"								"
" Optional parameters:						"
"	key=		Header key(s) to range (default=all)	"
" 								"
" Note: Gives partial results if interrupted			"
" 								"
" Output is: 							"
" number of traces 						"
" keyword min max (first - last) 				"
" 								"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace surange {


/* Credits:
 *	Geocon: Garry Perratt (output one header per line;
 *		option to specify headers to range;
 *		added first & last values where min<max)
 *	Based upon original by:
 *		SEP: Stew Levin
 *		CWP: Jack K. Cohen
 *
 * Note: the use of "signal" is inherited from BSD days and may
 *       break on some UNIXs.  It is dicy in that the responsibility
 *	 for program termination is lateraled back to the main.
 *
 */
/**************** end self doc ***********************************/


/* Prototypes */
void printrange(segy *tpmin, segy *tpmax, segy *tpfirst, segy *tplast);
static void closeinput(void);

segy tr, trmin, trmax, trfirst, trlast;

void* main_surange( void* args )
{
	int ntr;			/* number of traces		*/
	int nkeys=0;			/* number of keywords to range	*/
	Value val;			/* value of current keyword	*/
	Value valmin;			/* smallest seen so far		*/
	Value valmax;			/* largest seen so far		*/
	cwp_String type;		/* data type of keyword		*/
	cwp_String key[SU_NKEYS];	/* array of keywords		*/


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_surange );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get "key" value */
	if ((nkeys=parObj.countparval("key"))!=0) {
		parObj.getparstringarray("key",key);
	}

        parObj.checkpars();

	/* Zero out values of trmin and trmax */
	memset((void *) &trmin, 0, sizeof(segy));
	memset( (void *) &trmax, 0, sizeof(segy));

	/* Set up closing commands */
	signal(SIGINT, (void (*) (int)) closeinput);
	signal(SIGTERM, (void (*) (int)) closeinput);

	/* Do first trace outside loop to initialize mins and maxs */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	{	register int i;
		if (nkeys==0) {
   			for (i = 0; i < SU_NKEYS; ++i) {
 				gethval(&tr, i, &val);
 				puthval(&trmin, i, &val);
 				puthval(&trmax, i, &val);
 				puthval(&trfirst, i, &val);
			}
		} else	{
			register int j;
			for (i=0;i<nkeys;i++) {
				j = getindex(key[i]);
 				gethval(&tr, j, &val);
 				puthval(&trmin, j, &val);
 				puthval(&trmax, j, &val);
 				puthval(&trfirst, j, &val);
			}
		}
	}

	ntr = 1;
	while (cs2su->getTrace(&tr)) {
		register int i;
		if (nkeys==0) {
	       		for (i = 0; i < SU_NKEYS; ++i) {
				type = hdtype(getkey(i));
				gethval(&tr, i, &val);
				gethval(&trmin, i, &valmin);
				gethval(&trmax, i, &valmax);
				if (valcmp(type, val, valmin) < 0)
					puthval(&trmin, i, &val);
				if (valcmp(type, val, valmax) > 0)
					puthval(&trmax, i, &val);
 				puthval(&trlast, i, &val);
			}
		} else	{
			register int j;
			for (i=0;i<nkeys;i++) {
				type = hdtype(key[i]);
				j = getindex(key[i]);
				gethval(&tr, j, &val);
				gethval(&trmin, j, &valmin);
				gethval(&trmax, j, &valmax);
				if (valcmp(type, val, valmin) < 0)
					puthval(&trmin, j, &val);
				if (valcmp(type, val, valmax) > 0)
					puthval(&trmax, j, &val);
 				puthval(&trlast, j, &val);

			}
		}
		++ntr;
	}

	printf("%d traces:\n",ntr);
	printrange(&trmin, &trmax, &trfirst, &trlast);


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



/* printrange - print non-zero header values ranges	*/
void printrange(segy *tpmin, segy *tpmax, segy *tpfirst, segy *tplast)
{
	register int i = 0;
	Value valmin, valmax, valfirst, vallast;
	double dvalmin, dvalmax, dvalfirst, dvallast;
	cwp_String key;
	cwp_String type;
	int kmin = 0, kmax=SU_NKEYS;

	for (i = kmin; i < kmax; ++i) {
		key = getkey(i);
		type = hdtype(key);
		gethval(tpmin, i, &valmin);
		gethval(tpmax, i, &valmax);
		gethval(tpfirst, i, &valfirst);
		gethval(tplast, i, &vallast);
		dvalmin = vtod(type, valmin);
		dvalmax = vtod(type, valmax);
		dvalfirst = vtod(type, valfirst);
		dvallast = vtod(type, vallast);
		if (dvalmin || dvalmax) {
			if (dvalmin < dvalmax) {
				printf("%-8s ", key);
				printfval(type, valmin);
				printf(" ");
				printfval(type, valmax);
				printf(" (");
				printfval(type, valfirst);
				printf(" - ");
				printfval(type, vallast);
				printf(")");
			} else {
				printf("%-8s ", key);
				printfval(type, valmin);
			}
			putchar('\n');
		}
	}
	return;
}


static void closeinput(void) /* for graceful interrupt termination */
{
	/* Close stdin and open /dev/null in its place.  Now we are reading */
	/* from an empty file and the loops terminate in a normal fashion.  */

	efreopen("/dev/null", "r", stdin);
}

} // END namespace
