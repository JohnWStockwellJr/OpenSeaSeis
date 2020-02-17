/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUNHMOSPIKE: $Revision: 1.4 $ ; $Date: 2011/11/12 00:40:42 $	*/

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
std::string sdoc_sunhmospike =
"									"
" SUNHMOSPIKE - generates SPIKE test data set with a choice of several   "
"   Non-Hyperbolic MOveouts						"
"									"
"   sunhmospike [optional parameters] > out_data_file  			"
"									"
" Optional parameters:							"
"	nt=300	  number of time samples				"
"	ntr=20	  number of traces					"
"	dt=0.001	time sample rate in seconds			"
"	offref=2000	reference offset				"
"									"
"	gopt=		1 = parabolic transform model			"
"			2 = Foster/Mosher pseudo hyperbolic option model"
"			3 = linear tau-p model				"
"	depthref=400	reference depth used when gopt=2		"
"	offinc=100	offset increment				"
"	nspk=4	  number of events					"
"									"
"	p1 = 0		event moveout for event #1 in ms on reference offset"
"	t1 = 100	intercept time ms event #1			"
"	a1 = 1.0	amplitude for event #1				"
"									"
"	p2 = 200	event moveout for event #2 in ms on reference offset"
"	t2 = 100	intercept time ms for spike #2			"
"	a2 = 1.0	amplitude for event #2				"
"									"
"	p3 = 0;		event moveout for event #3 in ms on reference offset"
"	t3 = 200	intercept time for spike #3			"
"	a3 = 1.0	amplitude for event #3				"
"	p4 = 120	 event moveout for event #4 in ms on reference offset"
"									"
"	t4 = 200	intercept time s for spike #4			"
"	a4 = 1.0	amplitude for event #4				"
"									"
"	cdp = 1	 output cdp number					"
" 									"
" Notes:								"
" Creates a common cdp su data file with up to four spike events	"
" for impulse response studies for suradon, and sutifowler		"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sunhmospike {


/* Credits:
 *	CWP: Shuki Ronen, Chris Liner, 
 *      Modified: CWP   by John Anderson, April, 1994, to have
 *       appropriate trace header words and default values 
 *       for SUTIFOWLER tests
 *
 */
/**************** end self doc ***********************************/

float gofx(int igopt, float offset, float intercept_off,float refdepth);
segy tr;

void* main_sunhmospike( void* args )
{
	int ntr,nt,itr,igopt,j,cdp;
	int nspk,it1,it2,it3,it4;
	float dt,off,offinc;
	float t1,t2,t3,t4,p1,p2,p3,p4,a1,a2,a3,a4;
	float refdepth,offref;

	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sunhmospike );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	if(!parObj.getparint("gopt", &igopt)) igopt = 1;
	if(!parObj.getparint("ntr",&ntr)) ntr = 20;
	if(!parObj.getparint("nt", &nt)) nt = 300;
	if(!parObj.getparint("cdp",&cdp)) cdp=1;	

	dt = 0.001;	parObj.getparfloat("dt", &dt);		
	offref = 2000.;	parObj.getparfloat("offref", &offref);
	offinc = 100.;	parObj.getparfloat("offinc", &offinc);
	refdepth = 400.;parObj.getparfloat("refdepth", &refdepth);		
	nspk = 4;	parObj.getparint("nspk", &nspk);

	p1 = 0.;	parObj.getparfloat("p1", &p1);
	t1 = 100;	parObj.getparfloat("t1", &t1);
	a1 = 1.0;	parObj.getparfloat("a1", &a1);
	p2 = 200;	parObj.getparfloat("p2", &p2);
	t2 = 100;	parObj.getparfloat("t2", &t2);
	a2 = 1.0;	parObj.getparfloat("a2", &a2);
	p3 = 0.;	parObj.getparfloat("p3", &p3);
	t3 = 200;	parObj.getparfloat("t3", &t3);
	a3 = 1.0;	parObj.getparfloat("a3", &a3);
	p4 = 100;	parObj.getparfloat("p4", &p4);
	t4 = 200;	parObj.getparfloat("t4", &t4);
	a4 = 1.0;	parObj.getparfloat("a4", &a4);

	t1 = t1/(1000.*dt);
	t2 = t2/(1000.*dt);
	t3 = t3/(1000.*dt);
	t4 = t4/(1000.*dt);
	p1 = p1/(gofx(igopt,offref,0.,refdepth)*1000*dt);
	p2 = p2/(gofx(igopt,offref,0.,refdepth)*1000*dt);
	p3 = p3/(gofx(igopt,offref,0.,refdepth)*1000*dt);
	p4 = p4/(gofx(igopt,offref,0.,refdepth)*1000*dt);

	for( itr=0;itr<ntr; itr++) {
		for(j=0;j<nt;j++) tr.data[j]=0.;

		off = (itr+1)*offinc;
		it1 = t1+p1*gofx(igopt,off,0.,refdepth);
		it2 = t2+p2*gofx(igopt,off,0.,refdepth);
		it3 = t3+p3*gofx(igopt,off,0.,refdepth);
		it4 = t4+p4*gofx(igopt,off,0.,refdepth);

/*		fprintf(stderr,"itr=%d %d %d %d %d\n",itr,it1,it2,it3,it4); */

		if (nspk > 0 && it1>0 && it1<nt) tr.data[it1] += a1; 
		if (nspk > 1 && it2>0 && it2<nt) tr.data[it2] += a2;
		if (nspk > 2 && it3>0 && it3<nt) tr.data[it3] += a3;
		if (nspk > 3 && it4>0 && it4<nt) tr.data[it4] += a4;

		tr.tracl = itr + 1;
		tr.cdp = cdp;
		tr.offset = off;
		tr.ns = nt;
		tr.dt = dt*1000000;
		tr.f1=0;
		tr.d1=dt;
		tr.f2=offinc;
		tr.d2=offinc;

		su2cs->putTrace(&tr);
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
float gofx(int igopt, float offset, float intercept_off,float refdepth)
/*******************************************************************
return g(x) for various options
******************************************************************
Function parameters:

int igopt		1 = parabolic transform
			2 = Foster/Mosher pseudo hyperbolic option
			3 = linear tau-p
			4 = linear tau-p using absolute value of
				offset
float offset		offset in m
float intercept_off	offset corresponding to intercept time
float refdepth		reference depth in m for igopt=2
*******************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*******************************************************************/
{
	offset=offset-intercept_off;
	if(igopt==1) {
		return(offset*offset);
	}
	if(igopt==2) {
		return( sqrt(refdepth*refdepth + offset*offset) );
	}
	if(igopt==3) {
		return(offset);
	}
	if(igopt==4) {
		return(fabs(offset));
	}
	return(offset);
}

} // END namespace
