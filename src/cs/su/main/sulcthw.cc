/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SULCTHW: $Revision: 1.4 $ ; $Date: 2011/11/16 22:10:29 $	*/

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

#include "su.h"

/*********************** self documentation *****************************/
std::string sdoc_sulcthw =
" SULCTHW - Linear Coordinate Transformation of Header Words		"
"									"
"   sulcthw <infile >outfile						"
"									"
" xt=0.0	Translation of X					"
" yt=0.0	Translation of Y					"
" zt=0.0	Translation of Z					"
" xr=0.0	Rotation around X in degrees	 			"
" yr=0.0	Rotation aroun Y  in degrees	 			"
" zr=0.0	Rotation around Z in degrees 				"
"									"
" Notes:								"
" Translation:							"
" x = x'+ xt;y = y'+ yt;z = z' + zt;					"
"									"
" Rotations:					  			"
" Around Z axis							"
" X = x*cos(zr)+y*sin(zr);			  			"
" Y = y*cos(zr)-x*sin(zr);			  			"
" Around Y axis							"
" Z = z*cos(yr)+x*sin(yr);			  			"
" X = x*cos(yr)-z*sin(yr);			  			"
" Around X axis							"
" Y = y*cos(xr)+z*sin(xr);			  			"
" Z = Z*cos(xr)-y*sin(xr);			  			"
"									"
" Header words triplets that are transformed				"
" sx,sy,selev								"
" gx,gy,gelev								"
"									"
" The header words restored as double precision numbers using SEG-Y	"
" convention (with coordinate scalers scalco and scalel).		"
"									"
" After transformation they are converted back to			"
" short and stored, no.				 			"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sulcthw {


/*
 *  Credits: Potash Corporation of Saskatchewan: Balasz Nemeth   c. 2008
 *
 */

/**************** end self doc ********************************/
#define PP180 0.017453292

/* Segy data */
segy tr;				/* SEGY trace */

/* type defined to store coordinates */
typedef struct {
		double x;
		double y;
		double z;
} ctrp;

/* define transformation function */
ctrp transf(ctrp p,ctrp t,ctrp r);

void* main_sulcthw( void* args )
{

	/* declarations */
	ctrp p;
	ctrp pt;
	ctrp t;
	ctrp r;
	float tmp;

	/* Initargs */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sulcthw );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */

		
	
	/* Get parameters */
	if (!parObj.getparfloat("xt", &tmp)) tmp = 0.0;
	t.x = (double)tmp;
	if (!parObj.getparfloat("yt", &tmp)) tmp = 0.0;
	t.y = (double)tmp;
	if (!parObj.getparfloat("zt", &tmp)) tmp = 0.0;
	t.z = (double)tmp;
	if (!parObj.getparfloat("xr", &tmp)) tmp = 0.0;
	r.x = (double)tmp;
	if (!parObj.getparfloat("yr", &tmp)) tmp = 0.0;
	r.y = (double)tmp;
	if (!parObj.getparfloat("zr", &tmp)) tmp = 0.0;
	r.z = (double)tmp;
	
        parObj.checkpars();
	
	/* from degree to rad */
	r.x *= PP180;
	r.y *= PP180;
	r.z *= PP180;
	
	/* get info from first trace */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't get first trace");

	do {
		/* honor scalco and scalel */
		p.x = (double)tr.sx*pow(10.0,(double)tr.scalco);
		p.y = (double)tr.sy*pow(10.0,(double)tr.scalco);
		p.z = (double)tr.selev*pow(10.0,(double)tr.scalel);

		/* apply transformation */
		pt = transf(p,t,r);
		
		/* write new values to header words sx, sy, selev */
		tr.sx	=  NINT(pt.x/pow(10.0,(double)tr.scalco));
		tr.sy	=  NINT(pt.y/pow(10.0,(double)tr.scalco));
		tr.selev =  NINT(pt.z/pow(10.0,(double)tr.scalel));
		
		/* again, honor scalco */
		p.x = (double)tr.gx*pow(10.0,(double)tr.scalco);
		p.y = (double)tr.gy*pow(10.0,(double)tr.scalco);
		p.z = (double)tr.gelev*pow(10.0,(double)tr.scalel);
		
		/* transformation */
		pt = transf(p,t,r);
		
		/* write new values of gx,gy, gelev */
		tr.gx	=  NINT(pt.x/pow(10.0,(double)tr.scalco));
		tr.gy	=  NINT(pt.y/pow(10.0,(double)tr.scalco));
		tr.gelev =  NINT(pt.z/pow(10.0,(double)tr.scalel));
				
		/* output altered trace */
		su2cs->putTrace(&tr);

	} while(cs2su->getTrace(&tr));
	
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

ctrp transf(ctrp p,ctrp t,ctrp r)
/* linear transformation */
{

	ctrp pt;

	/* Translate */
	p.x +=t.x;
	p.y +=t.y;
	p.z +=t.z;
	
	/* Rotate */
	/* Z */
	pt.x = p.x*cos(r.z)+p.y*sin(r.z);
	pt.y = p.y*cos(r.z)-p.x*sin(r.z);
	pt.z = p.z;
	
	/* Y */
	pt.z = pt.z*cos(r.y)+pt.x*sin(r.y);
	pt.x = pt.x*cos(r.y)-pt.z*sin(r.y);
	
	/* X */
	pt.y = pt.y*cos(r.x)+pt.z*sin(r.x);
	pt.z = pt.z*cos(r.x)-pt.y*sin(r.x);
	
	return(pt);

}

} // END namespace
