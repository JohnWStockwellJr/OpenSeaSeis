/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSHW: $Revision: 1.25 $ ; $Date: 2011/11/16 22:10:29 $		*/

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

/*********************** self documentation ******************************/
std::string sdoc_sushw =
" 									"
" SUSHW - Set one or more Header Words using trace number, mod and	"
"	 integer divide to compute the header word values or input	"
"	 the header word values from a file				"
" 									"
" ... compute header fields						"
"   sushw <stdin >stdout key=cdp,.. a=0,..  b=0,.. c=0,.. d=0,.. j=..,.."
" 									"
" ... or read headers from a binary file				"
"   sushw <stdin > stdout  key=key1,..    infile=binary_file		"
" 									"
" 									"
" Required Parameters for setting headers from infile:			"
" key=key1,key2 ... is the list of header fields as they appear in infile"
" infile= 	binary file of values for field specified by		"
" 		key1,key2,...						"
" 									"
" Optional parameters ():						"
" key=cdp,...			header key word(s) to set 		"
" a=0,...			value(s) on first trace			"
" b=0,...			increment(s) within group		"
" c=0,...			group increment(s)	 		"
" d=0,...			trace number shift(s)			"
" j=ULONG_MAX,ULONG_MAX,...	number of elements in group		"
" 									"
" Notes:								"
" Fields that are getparred must have the same number of entries as key	"
" words being set. Any field that is not getparred is set to the default"
" value(s) above. Explicitly setting j=0 will set j to ULONG_MAX.	"
" 									"
" The value of each header word key is computed using the formula:	"
" 	i = itr + d							"
" 	val(key) = a + b * (i % j) + c * (int(i / j))			"
" where itr is the trace number (first trace has itr=0, NOT 1)		"
" 									"
" Examples:								"
" 1. set every dt field to 4ms						"
" 	sushw <indata key=dt a=4000 |...				"
" 2. set the sx field of the first 32 traces to 6400, the second 32 traces"
"    to 6300, decrementing by -100 for each 32 trace groups		"
"   ...| sushw key=sx a=6400 c=-100 j=32 |...				"
" 3. set the offset fields of each group of 32 traces to 200,400,...,6400"
"   ...| sushw key=offset a=200 b=200 j=32 |...				"
" 4. perform operations 1., 2., and 3. in one call			"
"  ..| sushw key=dt,sx,offset a=4000,6400,200 b=0,0,200 c=0,-100,0 j=0,32,32 |"
" 									"
" In this example, we set every dt field to 4ms.  Then we set the first	"
" 32 shotpoint fields to 6400, the second 32 shotpoint fields to 6300 and"
" so forth.  Next we set each group of 32 offset fields to 200, 400, ...,"
" 6400.									"
" 									"
" Example of a typical processing sequence using suchw:			"
"  sushw <indata key=dt a=4000 |					"
"  sushw key=sx a=6400 c=-100 j=32 |					"
"  sushw key=offset a=200 b=200 j=32 |			     		"
"  suchw key1=gx key2=offset key3=sx b=1 c=1 |		     		"
"  suchw key1=cdp key2=gx key3=sx b=1 c=1 d=2 >outdata	     		"
" 									"
" Again, it is possible to eliminate the multiple calls to both sushw and"
" sushw, as in Example 4.						"
" 									"
" Reading header values from a binary file:				"
" If the parameter infile=binary_file is set, then the values that are to"
" be set for the fields specified by key=key1,key2,... are read from that"
" file. The values are read sequentially from the file and assigned trace"
" by trace to the input SU data. The infile consists of C (unformated)	"
" binary floats in the form of an array of size (nkeys)*(ntraces) where	"
" nkeys is the number of floats in the first (fast) dimension and ntraces"
" is the number of traces.						"
" 									"
" Comment: 								"
" Users wishing to edit one or more header fields (as in geometry setting)"
" may do this via the following sequence:				"
"     sugethw < sudata output=geom  key=key1,key2 ... > hdrfile 	"
" Now edit the ASCII file hdrfile with any editor, setting the fields	"
" appropriately. Convert hdrfile to a binary format via:		"
"     a2b < hdrfile n1=nfields > binary_file				"
" Then set the header fields via:					"
"     sushw < sudata infile=binary_file key=key1,key2,... > sudata.edited"
" 									"
" Caveat: 								"
" If the (number of traces)*(number of key words) exceeds the number of	"
" values in the infile then the user may still set a single header field"
" on the remaining traces via the parameters key=keyword a,b,c,d,j.	"
"  									"
" Example:								"
"    sushw < sudata=key1,key2 ... infile=binary_file [Optional Parameters]"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sushw {


/* Credits:
 *	SEP: Einar Kajartansson
 *	CWP: Jack K. Cohen
 *      CWP: John Stockwell, added multiple fields and infile= options
 *
 * Caveat:
 *	All constants are cast to doubles.
 */
/**************** end self doc ****************************************/


segy tr;

/* Prototypes */
double mod(double x, double y);
void setval(cwp_String type, Value *valp, double a, double b,
		double c, double i, double j);


void* main_sushw( void* args )
{
	cwp_String key[SU_NKEYS];  /* array of keywords			*/
	cwp_String type[SU_NKEYS]; /* array of keywords			*/
	int index[SU_NKEYS];	/* name of type	of getparred key	*/

	int ikey;		/* key counter 				*/
	int nkeys;		/* number of header fields set		*/
	int count=0;		/* number of header fields from file	*/
	double i;		/* parameters for computing fields	*/
	int itr = 0;		/* trace counter 			*/
	Value val;		/* value of key field 			*/

	char *infile="";	/* name of input file of header values	*/
	FILE *infp=NULL;	/* pointer to input file		*/
	cwp_Bool from_file=cwp_false; /* is the data from infile?	*/

	float *afile=NULL;	/* array of "a" values from file	*/
	double *a=NULL;		/* array of "a" values			*/
	double *b=NULL;		/* array of "b" values			*/
	double *c=NULL;		/* array of "c" values			*/
	double *d=NULL;		/* array of "d" values			*/
	double *j=NULL;		/* array of "j" values			*/
	int n;			/* number of a,b,c,d,j values		*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sushw );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get "key" values */
	if ((nkeys=parObj.countparval("key"))!=0) {
		parObj.getparstringarray("key",key);

	} else {
		key[0]="cdp";
	}

	/* get types and indexes corresponding to the keys */
	for (ikey=0; ikey<nkeys; ++ikey) {
		type[ikey]=hdtype(key[ikey]);
		index[ikey]=getindex(key[ikey]);
	}


	/* get name of infile */
	parObj.getparstring("infile",&infile);

	/* if infile is specified get specified keys from file */
	if (*infile!='\0') {

		/* open infile */
		if((infp=efopen(infile,"r"))==NULL)
			throw cseis_geolib::csException("cannot open infile=%s\n",infile);

		/* set from_file flag */
		from_file=cwp_true;
	}

	/* If not from file, parObj.getpar a,b,c,d,j */
	if (!from_file) { 
		/* get "a" values */
		if ((n=parObj.countparval("a"))!=0) { 
			if (n!=nkeys)
			throw cseis_geolib::csException("number of a values not equal to number of keys");

			a=ealloc1double(n);
			parObj.getpardouble("a",a);
		} else {
			a=ealloc1double(nkeys);
			for (ikey=0; ikey<nkeys; ++ikey) a[ikey]=0.;
		}
		
		/* get "b" values */
		if ((n=parObj.countparval("b"))!=0) { 
			if (n!=nkeys)
			throw cseis_geolib::csException("number of b values not equal to number of keys");

			b=ealloc1double(n);
			parObj.getpardouble("b",b);
		} else {
			b=ealloc1double(nkeys);
			for (ikey=0; ikey<nkeys; ++ikey) b[ikey]=0.;
		}
		
		/* get "c" values */
		if ((n=parObj.countparval("c"))!=0) { 
			if (n!=nkeys)
			throw cseis_geolib::csException("number of c values not equal to number of keys");

			c=ealloc1double(n);
			parObj.getpardouble("c",c);
		} else {
			c=ealloc1double(nkeys);
			for (ikey=0; ikey<nkeys; ++ikey) c[ikey]=0.;
		}

		/* get "d" values */
		if ((n=parObj.countparval("d"))!=0) { 
			if (n!=nkeys)
			throw cseis_geolib::csException("number of d values not equal to number of keys");

			d=ealloc1double(n);
			parObj.getpardouble("d",d);
		} else {
			d=ealloc1double(nkeys);
			for (ikey=0; ikey<nkeys; ++ikey) d[ikey]=0.;
		}

		/* get "j" values */
		if ((n=parObj.countparval("j"))!=0) { 
			if (n!=nkeys)
			throw cseis_geolib::csException("number of j values not equal to number of keys");

			j=ealloc1double(n);
			parObj.getpardouble("j",j);

			/* make sure that j!=0 */
			for (ikey=0; ikey<nkeys; ++ikey)
				if(j[ikey]==0) j[ikey]=ULONG_MAX;
		} else {
			j=ealloc1double(nkeys);
			for (ikey=0; ikey<nkeys; ++ikey) j[ikey]=ULONG_MAX;
		}
	} else { /* if reading from a file */
		/* allocate space for afile */
		afile=ealloc1float(nkeys);
	}

        parObj.checkpars();

	/* loop over traces */
	while (cs2su->getTrace(&tr)) {

		if (from_file) {
			/* use the "a" value from file to trace by trace */
			if (efread(afile,FSIZE,nkeys,infp)!=0) {
				for (ikey=0; ikey<nkeys; ++ikey) {
					double a_in;
					a_in=(double) afile[ikey];
					setval(type[ikey],&val,a_in,
							 0,0,0,ULONG_MAX);
					puthval(&tr,index[ikey],&val);
				++count;
				}
			}
		} else { /* use parObj.getparred values of a,b,c,d,j */
			for (ikey=0; ikey<nkeys; ++ikey) {
				i = (double) itr + d[ikey];
				
				setval(type[ikey],&val,a[ikey],b[ikey],
						c[ikey],i,j[ikey]);
				puthval(&tr,index[ikey],&val);
			}

		}

		++itr;
		su2cs->putTrace(&tr);
	}

	if (from_file) {
		efclose(infp);
		if (count < (int)(itr*nkeys) ) {
		   warn("itr=%d > count=%d %s",(int) itr*count,count);
		   warn("n traces=%d > data count =%d",(itr*nkeys),count);
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


void setval( cwp_String type, Value *valp, double a, double b,
		double c, double i, double j)
{
	switch (*type) {
	case 's':
		throw cseis_geolib::csException("can't set char header word");
	break;
	case 'h':
		valp->h = a + b * mod(i, j) + c * ((int) (i/j));
	break;
	case 'u':
		valp->u = a + b * mod(i, j) + c * ((int) (i/j));
	break;
	case 'l':
		valp->l = (long) (a + b * mod(i, j) + c * ((int) (i/j)));
	break;
	case 'v':
		valp->v = (unsigned long) (a + b * mod(i, j) + c * ((int) (i/j)));
	break;
	case 'i':
		valp->i = a + b * mod(i, j) + c * ((int) (i/j));
	break;
	case 'p':
		valp->p = a + b * mod(i, j) + c * ((int) (i/j));
	break;
	case 'f':
		valp->f = a + b * mod(i, j) + c * ((int) (i/j));
	break;
	case 'd':
		valp->d = a + b * mod(i, j) + c * ((int) (i/j));
	default:
		throw cseis_geolib::csException("unknown type %s", type);
	break;
	}
	return;
}


double mod(double x, double y)	/* As defined in Knuth, vol. 1	*/
{
	return y == 0.0 ? x : x - y * floor(x/y);
}

} // END namespace
