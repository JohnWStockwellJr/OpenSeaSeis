/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* Copyright (c) Colorado School of Mines, 2010.*/
/* All rights reserved.                       */

/* SUAHW: $Revision: 1.3 $ ; $Date: 2011/11/17 00:03:38 $	*/

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
std::string sdoc_supad =
"									"
" SUPAD - Pad zero traces						"
"									"
"  supad <stdin >stdout min= max= [optional parameters]			"
"									"
" Required parameters:							"
"  min=			trace key start					"
"  max=			trace key end					"
"									"
" Optional parameters:							"
"  key1=ep		panel key 					"
"  key2=tracf		trace key 					"
"  key3=trid		flag key					"
"  val3=2		value assigned to padded traces			"
"  d=1			trace key spacing				"
"									"
" Notes:								"
"  In contrast to most SU codes, supad recognizes panels, or ensembles.	"
"  If the input consists of several panels, each panel will be padded	"
"  individually.							"
"  key1 and key2 are the primary and secondary sort key of the data set."
"  The sort order of key1 does not matter at all.			"
"  The sort order of key2 must be monotonous - if key2 is descending,	"
"	supply a negative value for the spacing d.			"
"  Traces with a key2-value outside the min/max range will be lost. 	"
"  Traces with a key2-value that is not a multiple of the spacing from	"
"	the min-value (the max-value, if the spacing is negative) will	"
"	not be lost. Instead, they will shift the series of key2-values."
"  By default the dead trace flag will be raised for the padded traces.	"
"  This should make it easy to remove the zero traces later on, if need be."
"									"
" Examples:								"
"	suplane | supad min=1 max=40 key1=offset key2=tracr | ...	"
"	... appends eight empty traces.					"
"									"
"	suplane | supad min=1 max=32 key1=offset key2=tracr d=0.5 | ...	"
"	... inserts a zero trace after each trace (even though the	"
"	header tracr is integer and cannot properly store the floats)	"
"									"
"	suplane | supad min=1 max=32 | ...				"
"	... produces an error because the panel and trace key are all 0."
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace supad {


/* Credits:
 *	Florian Bleibinhaus, U Salzburg, Austria
 */
/**************** end self doc ***********************************/

segy tr, nulltr;

/* Prototypes */
void assgnval(cwp_String type1, Value *valp1, double dval1);

void* main_supad( void* args )
{
	cwp_String key1,key2,key3;	/* panel/trace/flag key		*/
	cwp_String type1,type2,type3;	/* type for panel/trace/flag key*/
	int index1,index2,index3;	/* indexes for key1/2/3		*/
	Value val1,val2,val3;		/* value of key1/2/3		*/
	double dval1,dval2,dval3;	/* value of key1/2/3		*/
	double c;			/* trace key spacing		*/
	double dmin,dmax,dx;		/* trace key start/end/spacing	*/
	double panelno,traceno;
	int nt;				/* number of samples per trace	*/
	int isgn;			/* sort order			*/
	int iflag;			/* internal flag:		*/
					/* -1	exit			*/
					/* 0	regular mode		*/
					/* 1	first time		*/
					/* 2	trace out of range	*/

	/* initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_supad );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* get parameters */
	if (!parObj.getparstring("key1", &key1)) key1 = "ep";
	if (!parObj.getparstring("key2", &key2)) key2 = "tracf";
	if (!parObj.getparstring("key3", &key3)) key3 = "trid";
	if (!parObj.getpardouble("val3", &dval3)) dval3 = 2.;
	if (!parObj.getpardouble("d", &dx)) dx = 1.;
	if (!parObj.getpardouble("min", &dmin)) throw cseis_geolib::csException("need lower panel boundary MIN");
	if (!parObj.getpardouble("max", &dmax)) throw cseis_geolib::csException("need upper panel boundary MAX");
        parObj.checkpars();

	/* check parameters */
	if (dx==0) throw cseis_geolib::csException("trace spacing d cannot be zero");
	if (dmax<dmin) throw cseis_geolib::csException("max needs to be greater than min");

	if (dx<0) {
		isgn = -1;
	} else {
		isgn = 1;
	}

	/* get types and index values */
	type1  = hdtype(key1);
	type2  = hdtype(key2);
	type3  = hdtype(key3);
	index1 = getindex(key1);
	index2 = getindex(key2);
	index3 = getindex(key3);

	/* loop over traces */
	iflag = 1;
	while (iflag>=0) {

		if (cs2su->getTrace(&tr)) {
			/* get header values */
			gethval(&tr, index1, &val1);
			gethval(&tr, index2, &val2);
			dval1 = vtod(type1, val1);
			dval2 = vtod(type2, val2);
			/* Initialize zero trace */
			nt = tr.ns;
			memset( (void *) nulltr.data, 0, nt*FSIZE);
			if ( iflag==1 ) {
				panelno = dval1;
				traceno = dmin - dx;
			}
			iflag = 0;
			if ( dval2<dmin || dval2>dmax ) iflag = 2;
/*	fprintf(stderr,"if=%d, dmin=%8.0f, dmax=%8.0f\n",iflag,dmin,dmax);*/
		} else {
			iflag = -1;	/* exit flag */
		}
/*	fprintf(stderr,"if=%d, dval1=%8.0f, dval2=%8.0f\n",iflag,dval1,dval2);*/

		/* if new panel or last trace --> finish the previous panel */
		if ( panelno!=dval1 || iflag==-1 ) {
/*	fprintf(stderr,"finish previous\n");*/
			for (c=traceno+dx; isgn*c<=isgn*dmax; c=c+dx) {
					assgnval(type2, &val2, c);
					puthval(&nulltr, index2, &val2);
					assgnval(type3, &val3, dval3);
					puthval(&nulltr, index3, &val3);
					su2cs->putTrace(&nulltr);
			}
			traceno = dmin - dx;	/* reset to pad present panel */
		}

		/* if trace within boundaries --> pad the present panel */
		if ( iflag==0 ) {
/*	fprintf(stderr,"pad present, trn=%5.0f,dval2=%5.0f\n",traceno,dval2);*/
			memcpy( (void *) &nulltr, (const void *) &tr, 240);
			for (c=traceno+dx; isgn*c<isgn*dval2; c=c+dx) {
					assgnval(type2, &val2, c);
					puthval(&nulltr, index2, &val2);
					assgnval(type3, &val3, dval3);
					puthval(&nulltr, index3, &val3);
					su2cs->putTrace(&nulltr);
			}
		}

		/* write the present trace and save header indices */
		if ( iflag==0 ) {
			su2cs->putTrace(&tr);
			panelno = dval1;
			traceno = dval2;
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


void assgnval(cwp_String type1, Value *valp1, double dval1)
{
		switch (*type1) {
		case 's':
			throw cseis_geolib::csException("can't change char header word");
		break;
		case 'h':
			valp1->h = (short) dval1;
		break;
		case 'u':
			valp1->u = (unsigned short) dval1;
		break;
		case 'l':
			valp1->l = (long) dval1;
		break;
		case 'v':
			valp1->v = (unsigned long) dval1;
		break;
		case 'i':
			valp1->i = (int) dval1;
		break;
		case 'p':
			valp1->p = (unsigned int) dval1;
		break;
		case 'f':
			valp1->f = (float) dval1;
		break;
		case 'd':
			valp1->d = (double) dval1;
		break;
		default:
			throw cseis_geolib::csException("unknown type %s", type1);
		break;
		}
}

} // END namespace
