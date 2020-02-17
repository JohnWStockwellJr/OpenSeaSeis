/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

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

/*********************** self documentation ******************************/
std::string sdoc_surelanan =
"									"
" SURELANAN - REsiduaL-moveout semblance ANalysis for ANisotropic media	"
"									"
" surelan refl= npicks=    [optional parameters]			"
" 									"
" Required parameters:							"
" reflector file: reflec =						"
" number of points in the reflector file =				"
" 									"
" Optional Parameters:							"
" nr1=51		number of r1-parameter samples   		"
" dr1=0.01              r1-parameter sampling interval			"
" fr1=-0.25             first value of r1-parameter			"
" nr2=51		number of r2-parameter samples   		"
" dr2=0.01              r2-parameter sampling interval			"
" fr2=-0.25             first value of r2-parameter			"
" dzratio=5             ratio of output to input depth sampling intervals"
" nsmooth=dzratio*2+1   length of semblance num and den smoothing window"
" verbose=0             =1 for diagnostic print on stderr		"
" method=linear		for linear interpolation of the interface       "
" 			=mono for monotonic cubic interpolation of interface"
" 			=akima for Akima's cubic interpolation of interface "
" 			=spline for cubic spline interpolation of interface "
"									"
" Note: 								"
" 1. This program is part of Debashish Sarkar's anisotropic model building"
" technique. 								"
" 2. Input migrated traces should be sorted by cdp - surelan outputs a 	"
"    group of semblance traces every time cdp changes.  Therefore, the  "
"    output will be useful only if cdp gathers are input.  		"
" 3. The residual-moveout semblance for cdp gathers is based		"
"	on z(h)*z(h) = z(0)*z(0) + r1*h^2 + r2*h^4/[h^2+z(0)^2] where z "
"	depth and h is the half-offset.   				"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace surelanan {

/**************** end self doc *******************************************/

/* 
 * Trace header fields accessed:  ns, d1, f1, offset, cdp
 * Trace header fields modified:  ns, d1, offset, cdp
 */

segy tr;

void* main_surelanan( void* args )
{
	float r1;
	float r2;
	int nr1;	/* number of r-parameter samples */
	int nr2;	/* number of r-parameter samples */
	float dr1;	/* r-parameter samples sampling interval */
	float dr2;	/* r-parameter samples sampling interval */
	float fr1;	/* first r-parameter samples */
	float fr2;	/* first r-parameter samples */
	int ir1;	/* r-parameter samples index */
	int ir2;	/* r-parameter samples index */
	int dzratio;	/* ratio of output to input sampling intervals */
	int nsmooth;	/* length in samples of num and den smoothing window */
	int nz;		/* number of time samples per input trace */
	float dz;	/* time sampling interval for input traces */
	float fz;	/* time of first sample input and output */
	int iz;		/* input depth sample index */
	int izi;	/* depth sample index used in linear interpolation */
	int is;
	float zi;	/* normalized time for linear interpolation */
	float frac;	/* fractional distance from sample in interpolation */
	int gottrace;	/* =1 if an input trace was read */
	int done;	/* =1 if done with everything */
	int verbose;	/* =1 for diagnostic print */
	int npicks;	/* number of picks made on the reflector */
	int ipicks;
	long cdp;	/* cdp from current input trace header */
	long cdpprev;	/* cdp from previous input trace header */
	float half_offset;	/* offset from input trace header */
	float roffs2;
	float zn;	/* time after residual-moveout */
	float nsum;	/* semblance numerator sum */
	float dsum;	/* semblance denominator sum */
	float temp;	/* temporary scalar */
	float Z;
	float Zmin;
	float *data;	/*array[nz] of input trace */
	float *sem;	/*array[nr2] of semblance */
	float *xint,*zint,(*zind)[4]; /*arrays storing picked interface */
	float *xcdp,*zcdp;
	float ***num;	/*array[nr1][nr2][nz] of semblance numerators */
	float ***den;	/*array[nr1][nr2][nz] of semblance denominators */
	float ***nnz;	/*array[nr1][nr2][nz] for counting non-zero samples */
	char *refl;	/*name of the file containing pts picked on the refl.*/
	char *method="linear";
	FILE *dataptr;  /*pointer to the file containing pts picked on the refl.*/

	/* hook up getpar */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_surelanan );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* get parameters from the first trace */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");
	nz = tr.ns;
	dz = tr.d1;
 	fz = tr.f1;
	cdp = tr.cdp;
	half_offset = tr.offset/2;

	/* get required parameters */
	parObj.getparstring("refl",&refl);
	parObj.getparint("npicks",&npicks);

	/* get optional parameters */
	if (!parObj.getparint("nr1",&nr1)) nr1 = 51;
	if (!parObj.getparint("nr2",&nr2)) nr2 = 51;
	if (!parObj.getparfloat("dr1",&dr1)) dr1 = 0.01;
	if (!parObj.getparfloat("dr2",&dr2)) dr2 = 0.01;
	if (!parObj.getparfloat("fr1",&fr1)) fr1 = -0.25;
	if (!parObj.getparfloat("fr2",&fr2)) fr2 = -0.25;
	if (!parObj.getparint("dzratio",&dzratio)) dzratio = 5;
	if (!parObj.getparint("nsmooth",&nsmooth)) nsmooth = dzratio*2+1;
	if (!parObj.getparint("verbose",&verbose)) verbose = 0;
	parObj.getparstring("method",&method);
	
        parObj.checkpars();
	dataptr=fopen(refl,"r");
	
	/* allocate memory */
	data = ealloc1float(nz);
	num  = ealloc3float(nsmooth,nr2,nr1);
	den  = ealloc3float(nsmooth,nr2,nr1);
	nnz  = ealloc3float(nsmooth,nr2,nr1);
	sem  = ealloc1float(nr2);
	zint = ealloc1float(npicks);
	xint = ealloc1float(npicks);
	xcdp = ealloc1float(1);
        zcdp = ealloc1float(1);


                	xcdp[0] = cdp;

		/* Input picked interface */
                for(ipicks=0;ipicks<npicks;ipicks++)
                        fscanf(dataptr,"%f %f\n", &zint[ipicks], &xint[ipicks]);

                /* if linear interpolation or only one input sample */
                if (method[0]=='l') {
                                intlin(npicks,xint,zint,zint[0],zint[npicks-1],
                                1,xcdp,zcdp);
                /* else, if monotonic interpolation */
                } else if (method[0]=='m') {
                                zind = (float (*)[4])ealloc1float(npicks*4);
                                cmonot(npicks,xint,zint,zind);
                                intcub(0,npicks,xint,zind,1,xcdp,zcdp);

                /* else, if Akima interpolation */
                } else if (method[0]=='a') {
                                zind = (float (*)[4])ealloc1float(npicks*4);
                                cakima(npicks,xint,zint,zind);
                                intcub(0,npicks,xint,zind,1,xcdp,zcdp);

                /* else, if cubic spline interpolation */
                } else if (method[0]=='s') {
                                zind = (float (*)[4])ealloc1float(npicks*4);
                                csplin(npicks,xint,zint,zind);
                                intcub(0,npicks,xint,zind,1,xcdp,zcdp);

                /* else, if unknown method specified */
                } else {
                        throw cseis_geolib::csException("%s is an unknown interpolation method!\n",method);
                }

		Z = zcdp[0];
		warn("%f \n",Z);

	/* zero accumulators */
	for (ir1=0; ir1<nr1; ++ir1) 
	for (ir2=0; ir2<nr2; ++ir2) {
		for (iz=0; iz<nsmooth; ++iz) {
			num[ir1][ir2][iz] = 0.0;
			den[ir1][ir2][iz] = 0.0;
			nnz[ir1][ir2][iz] = 0.0;
		}
	}

	/* initialize flags */
	gottrace = 1;
	done = 0;

	/* remember previous cdp */
	cdpprev = tr.cdp;

	/* loop over input traces */
	do {

		/* if got a trace */
		if (gottrace) {

			/* determine offset and cdp */
			half_offset = tr.offset/2;
			cdp = tr.cdp;

			/* get trace samples */
			memcpy((void *) data, (const void *) tr.data,
				nz*sizeof(float));
		}

		/* if cdp has changed or no more input traces */
		if (cdp!=cdpprev || !gottrace) {

			/* set output trace header fields */
			tr.offset = 0;
			tr.cdp = (int) cdpprev;
			tr.ns = nr2;
			tr.d1 = dr2;

			/* loop over r-parameters */
			for (ir1=0; ir1<nr1; ++ir1) {
				for (ir2=0; ir2<nr2; ++ir2) {
				/* compute semblance quotients */
					nsum = dsum = 0.0;
					for (is=0; is<nsmooth; ++is) {
						nsum += num[ir1][ir2][is]*
							num[ir1][ir2][is];
						dsum += nnz[ir1][ir2][is]*
							den[ir1][ir2][is];
					}
					sem[ir2] = (dsum!=0.0?nsum/dsum:0.0);
					/*warn("%f \n",sem[ir2]);*/
				}
				/* output semblances */
				memcpy((void *) tr.data,(const void *) sem,
					nr2*sizeof(float));
				su2cs->putTrace(&tr);

			}

			/* diagnostic print */
			if (verbose) 
				warn("tsemblance output for cdp=%d\n",cdpprev);

			/* if no more input traces, break input trace loop */
			if (!gottrace) break;

			/* remember previous cdp */
			cdpprev = cdp;
		}

		/* loop over residual moveouts */
		for (ir1=0,r1=fr1; ir1<nr1; ++ir1,r1+=dr1) 
		for (ir2=0,r2=fr2; ir2<nr2; ++ir2,r2+=dr2) {
 			
			/* compute hyperbolic and nonhyperbolic terms */
			roffs2 = r1*half_offset*half_offset+r2*pow(half_offset,4.0)/(pow(half_offset,2.0)+pow(Z,2.0));

			/* determine mute depth after rmo */

			Zmin=Z-((nsmooth-1)*dz/2);
			/* do rmo via linear interpolation and  */
			/* accumulate semblance numerator and denominator */
			for (iz=0,zn=Zmin; iz<nsmooth; ++iz,zn+=dz) {
				temp = zn*zn+roffs2;
				zi = (temp>fz*fz)?(sqrt(temp)-fz)/dz:0;
				izi = zi; /* izi is `int' and zi is `float' */
				if (izi<nz-1) {
					frac = zi-izi;
					temp = (1.0-frac)*data[izi]+
						frac*data[izi+1];
					if (temp!=0.0) {
						num[ir1][ir2][iz] += temp;
						den[ir1][ir2][iz] += temp*temp;
						nnz[ir1][ir2][iz] += 1.0;
					}
			    	}
			}
		}
	
		/* get next trace (if there is one) */
		if (!cs2su->getTrace(&tr)) gottrace = 0;

	} while (!done);

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
