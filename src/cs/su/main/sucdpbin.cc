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

/*********************** self documentation **********************/
std::string sdoc_sucdpbin =
" 									"
" SUCDPBIN - Compute CDP bin number					"
" 									"
" sucdpbin <stdin >stdout xline= yline= dcdp=				"
" 									"
" Required parameters:							"
" xline=		array of X defining the CDP line		"
" yline=		array of Y defining the CDP line		"
" dcdp=			distance between bin centers			"
"									"
" Optional parameters							"
" verbose=0		<>0 output informations				"
" cdpmin=1001		min cdp bin number				"
" distmax=dcdp		search radius					"
" 									"
" xline,yline defines the CDP line made of continuous straight lines. 	"
" If a smoother line is required, use unisam to interpolate.		"
" Bin centers are located at dcdp constant interval on this line. 	"
" Each trace will be numbered with the number of the closest bin. If no  "
" bin center is found within the search radius. cdp is set to 0		"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sucdpbin {


/* Credits:
 * 2009 Dominique Rousset - Mohamed Hamza 
 *      Université de Pau et des Pays de l'Adour (France)
 */

/**************** end self doc ***********************************/

segy tr;

void* main_sucdpbin( void* args )
{
	float *xline=NULL;	/* array of key mute curve values	*/
	float *yline=NULL;	/* ...		mute curve time values 	*/
	float dcdp;	/* distance between bin centers */
	float *xbin=NULL, *ybin=NULL;
	float sl,ord,x2,a,b,c,delta,dx;
	float xmp, ymp, dist, distmin, distmax;
	int nxline, nyline, nbin;
	int verbose;
	int ipoint;
	int cdpmin;
	float scale;
	float length;
	int maxbin,ibin,iseg ;

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sucdpbin );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get parameters */
	if (!(nxline = parObj.countparval("xline")))
		throw cseis_geolib::csException("must give xline= vector");
	if (!(nyline = parObj.countparval("yline")))
		throw cseis_geolib::csException("must give yline= vector");
	if (nxline != nyline)
		throw cseis_geolib::csException("lengths of xline, yline must be the same");
	if (!parObj.getparfloat("dcdp", &dcdp)) throw cseis_geolib::csException("must give dcdp");

	/* allocate space */
	xline = ealloc1float(nxline);	parObj.getparfloat("xline", xline);
	yline = ealloc1float(nyline);	parObj.getparfloat("yline", yline);
	for (ipoint=1; ipoint<nxline; ++ipoint)
			if (xline[ipoint]<=xline[ipoint-1])
				throw cseis_geolib::csException("xline values must increase monotonically");


	/* Get verbose parameter */
	if (!parObj.getparint("verbose", &verbose)) verbose = 0;

	if (verbose) {
	warn ("%d points on the line",nxline);
	ipoint=0 ;
	do  {
		warn ("point # %d x= %f y= %f",ipoint, xline[ipoint],yline[ipoint]);
		++ipoint;
	} while (ipoint < nxline);
	}
	if (!parObj.getparint("cdpmin", &cdpmin)) cdpmin = 1001;
	if (!parObj.getparfloat("distmax", &distmax)) distmax = dcdp;

        parObj.checkpars();

	/* maximum length of the line */

	length=0;
	for (ipoint=1; ipoint<nxline; ++ipoint)
	length=length+sqrt(pow((xline[ipoint]-xline[ipoint-1]),2.0)+pow((yline[ipoint]-yline[ipoint-1]),2.0));
	maxbin=length/dcdp+1;
	xbin=ealloc1float(maxbin+1); ybin=ealloc1float(maxbin+1);

	ibin=0; xbin[0]=xline[0]; ybin[0]=yline[0];
	iseg=0; sl=(yline[1]-yline[0])/(xline[1]-xline[0]);
	ord=(yline[0]*xline[1]-yline[1]*xline[0])/(xline[1]-xline[0]);

	dx=dcdp*sqrt(1/(1+sl*sl));

	/* loop over segments */
	do {
		if (verbose) warn("iseg=%d",iseg);

		/* 1 segment binning */	
		do {
		
			++ibin; 
		
			xbin[ibin]=xbin[ibin-1]+dx; 
			ybin[ibin]=sl*xbin[ibin]+ord;

			if (verbose)
			warn("ibin=%d x=%f y=%f ",ibin,xbin[ibin],ybin[ibin]);

		} while (xbin[ibin]<= xline[iseg+1]); ++iseg; --ibin;
	
		/* parameters of the next segment */
		sl=(yline[iseg+1]-yline[iseg])/(xline[iseg+1]-xline[iseg]);
		dx=dcdp*sqrt(1/(1+sl*sl));
		ord=(yline[iseg]*xline[iseg+1]
			-yline[iseg+1]*xline[iseg])/(xline[iseg+1]-xline[iseg]);

		if (verbose) warn("sl=%f ord=%f",sl,ord);

		a=1+sl*sl; 
		b=-2*xbin[ibin]+2*sl*(ord-ybin[ibin]); 
		c=xbin[ibin]*xbin[ibin]+pow(ord-ybin[ibin],2.0)-dcdp*dcdp;
		delta=b*b-4*a*c;
		x2=(-b+sqrt(delta))/2/a;
		++ibin;
		xbin[ibin]=x2;
		ybin[ibin]=sl*xbin[ibin]+ord;
		if (verbose) 
			warn("delta=%f x2=%f ibin=%d xbin=%f",
					delta,x2,ibin,xbin[ibin]);
		if (verbose) warn("ibin=%d x=%f y=%f ",
				ibin,xbin[ibin],ybin[ibin]);	
	

	} while (xbin[ibin] <= xline[nxline-1]);
	nbin=ibin;
	

	if (verbose) {
	warn ("length of the line :%f ",length);
	warn ("maximum number of bin:%d",maxbin);
	warn ("actual number of bin:%d",nbin);
	}  
	
	/* Get info from first trace */
	if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("can't read first trace");


	/* Loop over traces */
	do {
		if (tr.scalco < 0 ) 
			scale=1./abs(tr.scalco); 
		else if (tr.scalco > 0)
			scale=tr.scalco;
		else {
			warn ("scalco = 0 ; 1 assumed") ;
			scale=1;
		}	
		xmp=(tr.gx+tr.sx)/2.*scale; 
		ymp=(tr.gy+tr.sy)/2.*scale;
		tr.cdp=cdpmin; 
		distmin=sqrt(pow(xmp-xbin[0],2.0)+pow(ymp-ybin[0],2.0));

		for (ibin=1 ; ibin < nbin ; ++ibin) {
			dist=sqrt(pow(xmp-xbin[ibin],2.0)+pow(ymp-ybin[ibin],2.0));
		/* warn ("ibin= %d dist= %f distmin= %f xmp= %f ymp= %f scale=%f", ibin, dist,distmin,xmp,ymp,scale); */
			if (dist < distmin) {
				distmin=dist;
				tr.cdp=ibin+cdpmin;
			}
		} 
		if (distmin > distmax) tr.cdp=0;
		su2cs->putTrace(&tr);
	} while (cs2su->getTrace(&tr));
	
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
