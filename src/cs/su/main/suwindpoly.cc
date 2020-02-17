/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUPLANE: $Revision: 1.7 $ ; $Date: 2011/11/17 00:03:38 $	*/

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
std::string sdoc_suwindpoly =
" 									"
" SUWINDPOLY - WINDow data to extract traces on or within a respective	"
"	POLYgonal line or POLYgon with coordinates specified by header	"
"	keyword values 							"
" 									"
"  suwindpoly <stdin [Required parameters] [Optional params] file=outfile"
" 									"
" Required parameters:							"
" x=x1,x2,...	list of X coordinates for vertices			"
" y=y1,y2,...	list of Y coordinates for vertices			"
" file=file1,file2,..	output filename(s)				"
"									"
" Optional parameters							"
" xkey=fldr	X coordinate header key					"
" ykey=ep	Y coordinate header key					"
" pass=0 	polyline mode: pass traces near the polygonal line	"
"		=1 pass all traces interior to polygon			"
"		=2 pass all traces exterior to polygon			"
"									"
" Optional parameters used in polyline pass=0 mode only:		"
" The following need to be given if the unit increments in the X & Y	"
" directions are not equal.  For example, if fldr increments by 1 and	"
" ep increments by 4 to form 25 x 25 m bins specify dx=25.0 & dy=6.25.	"
" The output binning key will be converted to integers by the scaling	"
" with the smaller of the two values.					"
"									"
" dx=1.0	unit increment distance in X direction			"
" dy=1.0	unit increment distance in Y direction			"
" ilkey=tracl	key for resulting inline index in polyline mode		"
" xlkey=tracr	key for resulting xline index in polyline mode		"
" dw=1.0	distance in X-Y coordinate units of extracted line	"
"		to pass points to output.  Width of resulting line is	"
"		2*dw.  Ignored if polygon mode is specified.		"
" Notes:								"
" In polyline mode (pass=0), a single trace may be output multiple times"
" if it meets the acceptance criteria (distance from line segment < dw)	"
" for multiple line segments. However, the headers will be distinct	"
" and will associate the output trace with a line segment. This		"
" behavior facilitates creation of 3D supergathers from polyline	"
" output. Use susort after running in polyline mode.			"
"									"
" x=& y=lists should be repeated for as many polygons as needed when  "
" pass=1 or pass=2. 							"
"									"
" In polygon mode, the polygon closes itself from the last vertex to	"
" the first.								"
"									"
" Example:								"
"  suwindpoly <input.su x=10,20,50 y=0,30,60 dw=10 pass=0 file=out.su	"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suwindpoly {


/*
 * Credits:  Reginald H. Beardsley	rhb@acm.org
 *	    originally: suxarb.c adapted from the SLT/SU package.
 */

/**************** end self doc ********************************/

/* Function prototypes of routines used internally */
void changeval(cwp_String type, Value *val, int f);
int InPolygon( float X0 ,float Y0 ,float *X ,float *Y ,int n );

/* Structure used internally */
int PolylineTransform( 
	 float x0	/* X coordinate of input point */
	,float y0	/* Y coordinate of input point */
	,float *x	/* X coordinates of polyline vertices */
	,float *y	/* Y coordinates of polyline vertices */
	,float *s	/* arc length coordinate of polyline vertices */
,int n	   /* number of vertices in polyline */
	,float *il	/* arc length coordinate list for point */
	,float *xl	/* crossline coordinate list for point */
);

/* Global variable */
segy tr;

void* main_suwindpoly( void* args )
{

	int i;			/* counter				*/
	int j;			/* counter				*/

	cwp_String xkey="fldr";	/* x-coordinate keyword			*/
	cwp_String xtype;	/* ... its type				*/
	Value  xval;		/* ... its value			*/
	int    xindex;		/* ... its index			*/

	cwp_String ykey="ep";	/* y-coordinate keyword			*/
	cwp_String ytype;	/* ... its type				*/
	Value  yval;		/* ... its value			*/
	int    yindex;		/* ... its index			*/
	
	cwp_String ilkey="tracl";	/* in-line coordinate keyword	*/
	cwp_String iltype;		/* ... its type			*/
	Value  ilval;			/* ... its value		*/
	int    ilindex;			/* ... its index		*/

	cwp_String xlkey="tracr";	/* x-line coordinate keyword	*/
	cwp_String xltype;		/* ... its type			*/
	Value  xlval;			/* ... its value		*/
	int    xlindex;			/* ... its index		*/
	
	cwp_String *outfile=NULL;	/* outputfilename */

	float dx=1.0;	/* x sampling interval */
	float dy=1.0;	/* y sampling interval */

	float dw=1.0;	/* distance in X-Y coordinate units of	*/
			/* extracted line to pass points to output.*/
			/* Width of resulting line is 2*dw.	*/
			/* Ignored if polygon mode is specified.*/
	int nx=0; 		/* number of user-specified x-values	*/
	int ny=0;		/* number of user specifd y-values	*/
	
	int nFiles=0; 		/* number of output files specified */
	FILE **outfp=NULL;	/* filepointer to output file */

	int nPolygons=0;	/* number of polygons specified */
	int pass=0;		/* flag to specify the mode of the program */

	char errString[1024];	/* error message */

	float **xpoly=NULL;	/* array of x values defining polygon */
	float **ypoly=NULL;	/* array of y values defining polygon */
	float **spoly=NULL;	/* array of polygon diagonal distances */
	int    *npoly=NULL;	/* array of numbers of polygons */

	float x0=0.0;		/* first x value */
	float y0=0.0;		/* first y value */
	float *il=NULL;		/* in-line coordinates */
	float *xl=NULL;		/* x-line coordinates	*/
	int nl=0;
	int inside=0;


	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suwindpoly );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get parameters */
	parObj.getparstring("xkey",&xkey); 
	parObj.getparstring("ykey",&ykey); 

	parObj.getparstring("ilkey",&ilkey); 
	parObj.getparstring("xlkey",&xlkey); 

	if (!parObj.getparfloat("dw",&dw))		dw=1.0;
	if (!parObj.getparfloat("dx",&dx))		dx=1.0;
	if (!parObj.getparfloat("dy",&dy))		dy=1.0;

	if (!parObj.getparint("pass",&pass))		pass=0;
        parObj.checkpars();

	/* Count the number of x and y coordinates specified */
	nx=parObj.countparname("x");
	ny=parObj.countparname("y");

	/* Error trapping */
	if (pass==0 ){
		if (nx>1 || ny>1) {
			throw cseis_geolib::csException("only one set of x=& y=permitted if pass=0");

		} else if (nx!=1) {
			throw cseis_geolib::csException( "x values missing for polyline");

		} else if(ny!=1) {
			throw cseis_geolib::csException("y values missing for polyline");
		} else {
				nPolygons=1;
		}

	} else if (nx>ny){
		throw cseis_geolib::csException( "y values missing for polygon" );
	} else if (ny>nx){
		throw cseis_geolib::csException( "x values missing for polygon" );

	} else {
		nPolygons=nx;
	}

	/*----------------------*/
	/* get output filenames */
	/*----------------------*/
	nFiles=parObj.countparname( "file" );

	if( pass==1 &&  nFiles !=nPolygons ){
		throw cseis_geolib::csException( "One output filename must be given for each polygon" );

	} else if(pass==2 && nFiles !=1 ){
		throw cseis_geolib::csException( "Only one output filename may be given when pass=2" );

	} else if( pass > 2 || pass < 0 ){
		throw cseis_geolib::csException( "Invalid pass=option" );
	
	}

	if( !(outfile=calloc(nFiles ,sizeof(char*) )) ){
		sprintf(errString,
			"calloc(3c) failed %s %d" ,__FILE__ ,__LINE__ );
		throw cseis_geolib::csException(errString);

	}

	if( !(outfp=calloc(nFiles ,sizeof(FILE*) )) ){
		sprintf(errString 
			,"calloc(3c) failed %s %d" ,__FILE__ ,__LINE__ );
		throw cseis_geolib::csException(errString);

	}

	/* Get the names and a count of outfiles */
	for(i=0; i<nFiles;++i){
		parObj.getnparstring(i+1 ,"file" ,&outfile[i] );
	}

	/*---------------------------------------*/
	/* allocate space for polygon boundaries */
	/*---------------------------------------*/
	if( !(xpoly=calloc(nPolygons ,sizeof(float*) )) ){
		sprintf(errString,
			"calloc(3c) failed %s %d" ,__FILE__ ,__LINE__ );
		throw cseis_geolib::csException(errString );
	
	} else if( !(ypoly=calloc(nPolygons ,sizeof(float*) )) ){
		sprintf(errString,
			"calloc(3c) failed %s %d" ,__FILE__ ,__LINE__ );
		throw cseis_geolib::csException(errString);

	} else if( !(spoly=calloc(nPolygons ,sizeof(float*) )) ){
		sprintf(errString,
			"calloc(3c) failed %s %d" ,__FILE__ ,__LINE__ );
		throw cseis_geolib::csException(errString);

	} else if( !(npoly=calloc( nx ,sizeof(int) )) ){
		sprintf(errString,
			"calloc(3c) failed %s %d" ,__FILE__ ,__LINE__ );
		throw cseis_geolib::csException(errString);


	}
	
	/*---------------------------------------------*/
	/* allocate space for vertices of each polygon */
	/*---------------------------------------------*/
	for(i=0; i<nPolygons;++i){

		nx=parObj.countnparval( i+1 ,"x" ); 
		ny=parObj.countnparval( i+1 ,"y" ); 

		if( nx > ny ){
			sprintf(errString,"Y values missing from polygon %d" ,i+1 );
			throw cseis_geolib::csException(errString);

		} else if( ny > nx ){
			sprintf(errString,"X values missing from polygon %d" ,i+1 );
			throw cseis_geolib::csException(errString);

		} else if( !(xpoly[i]=calloc( nx+1 ,sizeof(float) )) ){
			sprintf(errString
	 ,"calloc(3c) failed %s %d" ,__FILE__ ,__LINE__ );
			throw cseis_geolib::csException(errString);

		} else if( !(ypoly[i]=calloc( nx+1 ,sizeof(float) )) ){
			sprintf(errString
	 ,"calloc(3c) failed %s %d" ,__FILE__ ,__LINE__ );
			throw cseis_geolib::csException(errString);

		} else if( !(spoly[i]=calloc( nx+1 ,sizeof(float) )) ){
			sprintf(errString
	 ,"calloc(3c) failed %s %d" ,__FILE__ ,__LINE__ );
			throw cseis_geolib::csException(errString);

		} else {

			/*-----------------------*/
			/* read polygon vertices */
			/*-----------------------*/

			npoly[i]=nx;
			parObj.getnparfloat( i+1 ,"x" ,xpoly[i] );
			parObj.getnparfloat( i+1 ,"y" ,ypoly[i] );

			/*-------------------------*/
			/* close polygon if needed */
			/*-------------------------*/

			if( pass !=0 

			   && (xpoly[i][nx-1] !=xpoly[i][0]
			    || ypoly[i][ny-1] !=ypoly[i][0])  ){

			   xpoly[i][nx]=xpoly[i][0];
			   ypoly[i][ny]=ypoly[i][0];
			   npoly[i]=nx + 1;

			} else if( !(il=calloc( nx ,sizeof(float) )) ){
			   sprintf(errString
	    ,"calloc(3c) failed %s %d" ,__FILE__ ,__LINE__ );
			   throw cseis_geolib::csException(errString);

			} else if( !(xl=calloc( nx ,sizeof(float) )) ){
			   sprintf(errString
	    ,"calloc(3c) failed %s %d" ,__FILE__ ,__LINE__ );
			   throw cseis_geolib::csException(errString);

			} else {

			   /*-------------*/
			   /* scale units */
			   /*-------------*/

			   for(j=0; j<npoly[i]; j++ ){

				xpoly[i][j] *=dx;
				ypoly[i][j] *=dy;

			   }

			}

			/*----------------------------------*/
			/* calculate arc length to vertices */
			/*----------------------------------*/

			for(j=1; j<npoly[i]; j++ ){

			   spoly[i][j]=spoly[i][j-1]
		+ sqrt( pow((xpoly[i][j]-xpoly[i][j-1]) ,2.0) 
			+ pow((ypoly[i][j]-ypoly[i][j-1]) ,2.0) );

			}
		}
	}

	/*---------------------------------------------------*/
	/* report the polyline or polygon vertices tostderr*/
	/*---------------------------------------------------*/

	for(i=0; i<nPolygons;++i){

		if( pass==1 ){
			fprintf(stderr,"Passing interior of ");
		} else if( pass==2 ){
			fprintf(stderr ,"Passing exterior of ");
		}
fprintf(stderr,"Polygon: %4d " ,i+1	 );
		if( nFiles > i ){
			fprintf(stderr,"Filename: %s " ,outfile[i] );
		}

		warn("");

		for(j=0; j<npoly[i]; j++ ){

fprintf(stderr,"%14s"	 ,""	   );
fprintf(stderr,"Vertex: %d " ,j+1	  );
			fprintf(stderr,"X: %f "	,xpoly[i][j] );
			fprintf(stderr,"Y: %f "	,ypoly[i][j] );
			fprintf(stderr,"S: %f "	,spoly[i][j] );

			fprintf(stderr,"\n" );
		}
		fprintf(stderr,"\n" );
	}

	/*-----------------------*/
	/* open the output files */
	/*-----------------------*/

	for(i=0; i<nFiles;++i){
		if( !(outfp[i]=fopen( outfile[i] ,"w" )) ){
			sprintf(errString,
				"fopen(3c) failed %s %d" ,__FILE__ ,__LINE__ );
			throw cseis_geolib::csException(errString);
		}

	}

	/*---------------------*/
	/* read in first trace */
	/*---------------------*/

	/* Get the first trace */
	if (!cs2su->getTrace(&tr) ){
		throw cseis_geolib::csException("can't get first trace");
	}

	/* Get types and indexes */
	xtype=hdtype(xkey);
	xindex=getindex(xkey);

	ytype=hdtype(ykey);
	yindex=getindex(ykey);

	iltype=hdtype(ilkey);
	ilindex=getindex(ilkey);

	xltype=hdtype(xlkey);
	xlindex=getindex(xlkey);
	
	/*----------------------------*/
	/* process all the input data */
	/*----------------------------*/
	do {
	
		gethval(&tr, xindex, &xval);
		x0=vtof(xtype, xval);

		gethval(&tr, yindex, &yval);
		y0=vtof(ytype, yval);

		if( pass==0 ){

			x0*=dx;
			y0*=dy;

			/*---------------------------------*/
			/* polyline transform & extraction */
			/*---------------------------------*/

			nl=PolylineTransform( x0 
		    ,y0 
		    ,xpoly[0] 
		    ,ypoly[0] 
		    ,spoly[0]
		    ,npoly[0]
		    ,il
		    ,xl
		   );

			for(j=0; j<nl; j++ ){

			   if( fabs(xl[j]) <=dw ){

				changeval( iltype ,&ilval ,(int) il[j] );
				puthval( &tr ,ilindex ,&ilval );

				changeval( xltype ,&xlval ,(int) xl[j] );
				puthval( &tr ,xlindex ,&xlval );
			   
				fputtr( outfp[0] ,&tr );
			   }

			}

		} else {

			/*--------------------*/
			/* polygon extraction */
			/*--------------------*/

			if( pass==1 ){

			   for(i=0; i<nPolygons;++i){

				inside=InPolygon( x0 
			    ,y0 
			    ,xpoly[i] 
			    ,ypoly[i] 
			    ,npoly[i] 
			   );

				if( inside ){
	  fputtr( outfp[i] ,&tr );
  
				}

			   }

			} else if( pass==2 ){

			   for(i=0; i<nPolygons;++i){

				inside=InPolygon( x0 
			    ,y0 
			    ,xpoly[i] 
			    ,ypoly[i] 
			    ,npoly[i] 
			   );

				if( inside ){
	  goto next;
				}

			   }

			   fputtr( outfp[0] ,&tr );

			}

		}

		next:;
			
	} while (fgettr(stdin ,&tr));


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

#include <stdio.h>
int InPolygon( float X0 ,float Y0 ,float* X ,float* Y ,int n )
/*--------------------------------------------------------------------*\
InPolygon() - determines if a point is within a polygon denoted by
 a polyline with identical first and last points.

	The basic algorithm is from "Computational Geometry in C", by
	Joseph O'Rourke but has been completely rewritten to make it 
	suitable for industrial use.

	The function returns integer values indicating the location of
	the input point relative to the bounding polygon:

	0 - point is exterior to the polygon
	1 - point is a vertex of the polygon
	2 - point lies on a line segment
	3 - point is strictly interior to the polygon

Reginald H. Beardsley				rhb@acm.org
\*--------------------------------------------------------------------*/
{

	int i;		/* vertex	   */
	int j;		/* adjacent vertex */

	int R=0;	/* number of right crossings */
	int L=0;	/* number of left crossings  */

	float x;

	for(i=0; i<n;++i){
		if( X[i]==X0 && Y[i]==Y0 ){
			return( 1 );
		}

		j=(i + n - 1) % n;

		if(   ((Y[i] > Y0) && (Y[j] <=Y0)) 
			 ||((Y[j] > Y0) && (Y[i] <=Y0)) ){
			
			x=((X[i]-X0) * (Y[j]-Y0) - (X[j]-X0) * (Y[i]-Y0))
			   / ((Y[j]-Y0) - (Y[i]-Y0));

			if( x > 0 ){
			   R++;
			}
		}

		if(   ((Y[i] > Y0) && (Y[j] <=Y0)) 
			 ||((Y[j] > Y0) && (Y[i] <=Y0)) ){
			
			x=((X[i]-X0) * (Y[j]-Y0) - (X[j]-X0) * (Y[i]-Y0))
			   / ((Y[j]-Y0) - (Y[i]-Y0));

			if( x < 0 ){
			   L++;
			}
		}

	}

	if( (R % 2) !=(L % 2) ){
		return( 2 ); 
	}

	if( (R % 2)==1 ){
		return( 3 );

	} else {
		return(0);

	}
}



#include <math.h>
#include <string.h>

int
PolylineTransform( 
	 float x0	/* X coordinate of input point */
	,float y0	/* Y coordinate of input point */
	,float* x	/* X coordinates of polyline vertices */
	,float* y	/* Y coordinates of polyline vertices */
	,float* s	/* arc length coordinate of polyline vertices */
	,int n		/* number of vertices in polyline */
	,float* il	/* arc length coordinate list for point */
	,float* xl	/* crossline coordinate list for point */
	)
/*====================================================================*\
PolylineTransform() calculates the position of the input point in
	arc length - arc normal coordinates for the polyline specified by
	the input vertex list.  It returns the number of segments for which
	the point lies on a normal to the segment.

	The coordinates are calculated relative to all segments for which 
	input point lies on a normal to the segment.

	The line direction metric starts at the first point.  The line normal
	metric is positive on the right side of the line.

Reginald H. Beardsley				rhb@acm.org
\*====================================================================*/
{

	float a;   /* intercept of line */
	float b;   /* slope of line	*/

	float x_;  /* X coordinate of intersection of normal & segment */
	float y_;  /* Y coordinate of intersection of normal & segment */

	float m_;  /* tangent of angle between m1 & m2			*/
	float m1;  /* vector from first point on segment to intersection */
	float m2;  /* vector from first point on segment to point	 */

	float dx;  /* distance of point from line segement along normal  */

	int i;
	int nl=0;

	memset( il ,0 ,sizeof(float)*n );
	memset( xl ,0 ,sizeof(float)*n );

	for(i=1; i<n;++i){

		if( x[i] !=x[i-1] ){

			/*---------------------------*/
			/* non-vertical line segment */
			/*---------------------------*/
			b=(y[i] - y[i-1]) / (x[i] - x[i-1] );
			a=y[i-1] - x[i-1] * b;
			x_=(x0 + (y0 - a) * b) / (1.0 + b*b);
			y_=a + b * x_;

			m1=(y_ - y[i-1]) / (x_ - x[i-1]);
			m2=(y0 - y[i-1]) / (x0 - x[i-1]);
	
			m_=(m1 - m2) / (1.0 + m1*m2);

			if( (x[i] >=x_ && x_ >=x[i-1]) 
				|| (x[i] <=x_ && x_ <=x[i-1]) ){

			   il[nl]=s[i-1] + sqrt( (y_-y[i-1])*(y_-y[i-1]) 
			    + (x_-x[i-1])*(x_-x[i-1]) );

			   dx=sqrt( (x0-x_)*(x0-x_) + (y0-y_)*(y0-y_) );

			   if( m_ >=0.0 ){
				xl[nl]=dx;

			   } else {
				xl[nl]=dx * -1.0;

			   }

			   nl++;
			}

		} else {

			/*-----------------------*/
			/* vertical line segment */
			/*-----------------------*/
			if ((y[i] >=y0 && y0 >=y[i-1]) 
				|| (y[i] <=y0 && y0 <=y[i-1]) ){
			
			   x_=x[i];
			   y_=y0;

			   il[nl]=s[i-1] + fabs( y_ - y[i-1] );

			   dx=fabs( x0 - x_);

			   if( y0 >=y[i-1] && x0 > x_ ){
				xl[nl]=dx * -1.0;

			   } else {
				xl[nl]=dx;

			   }
			   nl++;
			}
		}
	}

	return(nl);

}


void changeval(cwp_String type, Value *val, int f) {

	switch (*type) {
		  case 's':
				throw cseis_geolib::csException("can't change char header word");
		  break;
		  case 'h':
				val->h=f;
		  break;
		  case 'u':
				val->u=f;
		  break;
		  case 'l':
				val->l=f;
		  break;
		  case 'v':
				val->v=f;
		  break;
		  case 'i':
				val->i=f;
		  break;
		  case 'p':
				val->p=f;
		  break;
		  case 'f':
				val->f=f;
		  break;
		  case 'd':
				val->d=f;
		  break;
		  default:
				throw cseis_geolib::csException("unknown type %s", type);
		  break;
		  }
}

} // END namespace
