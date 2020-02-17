/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUAZIMUTH: $Revision: 1.20 $ ; $Date: 2011/11/16 22:10:29 $        */

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
std::string sdoc_suazimuth =
"                                                                        "
" SUAZIMUTH - compute trace AZIMUTH, offset, and midpoint coordinates    "
"             and set user-specified header fields to these values       "
"                                                                        "
"  suazimuth <stdin >stdout [optional parameters]                        "
"                                                                        "
" Required parameters:                                                   "
"     none                                                               "
"                                                                        "
" Optional parameters:                                                   "
"   key=otrav      header field to store computed azimuths in            "
"   scale=1.0      value(key) = scale * azimuth                          "
"   az=0           azimuth convention flag                               "
"                   0: 0-179.999 deg, reciprocity assumed                "
"                   1: 0-359.999 deg, points from receiver to source     "
"                  -1: 0-359.999 deg, points from source to receiver     "
"   sector=1.0     if set, defines output in sectors of size             "
"                  sector=degrees_per_sector, the default mode is        "
"                  the full range of angles specified by az              "
"                                                                        "
"   offset=0       if offset=1 then set offset header field              "
"   offkey=offset  header field to store computed offsets in             "
"                                                                        "
"   cmp=0          if cmp=1, then compute midpoint coordinates and       "
"                  set header fields for (cmpx, cmpy)                    "
"   mxkey=ep       header field to store computed cmpx in                "
"   mykey=cdp      header field to store computed cmpy in                "
"                                                                        "
" Notes:                                                                 "
"   All values are computed from the values in the coordinate fields     "
"   sx,sy (source) and gx,gy (receiver).                                 "
"   The output field \"otrav\" for the azimuth was chosen arbitrarily as "
"   an example of a little-used header field, however, the user may      "
"   choose any field that is convenient for his or her application.      "
"                                                                        "
"   Setting the sector=number_of_degrees_per_sector sets key field to    "
"   sector number rather than an angle in degrees.                       "
"                                                                        "
"   For az=0, azimuths are measured from the North, however, reciprocity "
"   is assumed, so azimuths go from 0 to 179.9999 degrees. If sector     "
"   option is set, then the range is from 0 to 180/sector.               "
"                                                                        "
"   For az=1, azimuths are measured from the North, with the assumption  "
"   that the direction vector points from the receiver to the source.    "
"   For az=-1, the direction vector points from the source to the        "
"   receiver. No reciprocity is assumed in these cases, so the angles go "
"   from 0 to 359.999 degrees.                                           "
"   If the sector option is set, then the range is from 0 to 360/sector. "
"                                                                        "
" Caveat:                                                                "
"   This program honors the value of scalco in scaling the values of     "
"   sx,sy,gx,gy. Type \"sukeyword scalco\" for more information.         "
"                                                                        "
"   Type \"sukeyword -o\" to see the keywords and descriptions of all    "
"   header fields.                                                       "
"                                                                        "
"   To plot midpoints, use: su3dchart                                    "
"                                                                        "
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suazimuth {


/* Credits:
 *  based on suchw, su3dchart
 *      CWP: John Stockwell and  UTulsa: Chris Liner, Oct. 1998
 *      UTulsa: Chris Liner added offset option, Feb. 2002
 *         cll: fixed offset option and added cmp option, May 2003
 *      RISSC: Nils Maercklin added key options for offset and 
 *             midpoints, and added azimuth direction option, Sep. 2006
 *
 *  Algorithms:
 *      offset = osign * sqrt( (gx-sx)*(gx-sx) + (gy-sy)*(gy-sy) )
 *               with osign = sgn( min((sx-gx),(sy-gy)) )
 *
 *      midpoint x  value  xm = (sx + gx)/2
 *      midpoint y  value  ym = (sy + gy)/2
 * 
 *  Azimuth will be defined as the angle, measured in degrees,
 *  turned from North, of a vector pointing to the source from the midpoint, 
 *  or from the midpoint to the source. Azimuths go from 0-179.000 degrees
 *  or from 0-180.0 degrees.
 *   
 *  value(key) = scale*[90.0 - (180.0/PI)*(atan((sy - ym)/(sx - xm))) ]
 *      or
 *  value(key) = scale*[180.0 - (180.0/PI)*(atan2((ym - sy),(xm - sx)) ]
 * 
 *  Trace header fields accessed: sx, sy, gx, gy, scalco. 
 *  Trace header fields modified: (user-specified keys)
 *
 */
/**************** end self doc ***********************************/

/* Prototype of functions used internally */
double computeAzimuth(double sx, double sy, 
    double mx, double my, double scale, int az);
void setval(cwp_String type_out, Value *val_out, double dval_out);

/* SEG-Y trace */
segy tr;

void* main_suazimuth( void* args )
{
    cwp_String key;        /* output key for computed azimuth */
    cwp_String type_out;   /* ... type for output key */
    int index_out;         /* ... index for output key */
    Value val_out;         /* ... value of output key */
    cwp_String okey="";       /* output key for computed offset */
    cwp_String otype_out="";  /* ...type for output key */
    int oindex_out=0;        /* ...index for output key */
    cwp_String mxkey="";      /* output key for computed cmp x */
    cwp_String mxtype_out=""; /* ...type for output key */
    int mxindex_out=0;       /* ...index for output key */
    cwp_String mykey="";      /* output key for computed cmp y */
    cwp_String mytype_out=""; /* ...type for output key */
    int myindex_out=0;       /* ...index for output key */

    double sx=0.0;         /* source x coordinate */
    double sy=0.0;         /* source y coordinate */
    double gx=0.0;         /* receiver x coordinate */
    double gy=0.0;         /* receiver y coordinate */
    double factor=0.0;     /* scaling factor */
    double scale=0.0;      /* value of scale */

    double mx, my;         /* x,y midpoint coordinates */
    int az=0;              /* azimuth convention flag */
    double sector=0;       /* azimuth by sectors */
    double azimuth;        /* computed azimuth value */
    int oflag=0;           /* offset flag */
    double offset;         /* computed offset value */
    int cmpflag=0;         /* cmp flag */
    double cmpx,cmpy;      /* cmp coords */
    
    /* Initialize */
    cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
    cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
    cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
    int argc = suArgs->argc;
    char **argv = suArgs->argv;
    cseis_su::csSUGetPars parObj;

    void* retPtr = NULL;  /*  Dummy pointer for return statement  */
    su2cs->setSUDoc( sdoc_suazimuth );
    if( su2cs->isDocRequestOnly() ) return retPtr;
    parObj.initargs(argc, argv);

    try {  /* Try-catch block encompassing the main function body */


    /* Get parameters */
    if (!parObj.getparstring("key",&key))         key = "otrav";
    if (!parObj.getparstring("offkey",&okey))   okey = "offset";
    if (!parObj.getparstring("mxkey",&mxkey))   mxkey = "ep";
    if (!parObj.getparstring("mykey",&mykey))   mykey = "cdp";
    if (!parObj.getpardouble("scale",&scale))     scale = 1.0;
    if (!parObj.getparint("az",&az))        az = 0;
    if (!parObj.getpardouble("sector",&sector))     sector = 1.0;
    if (!parObj.getparint("offset",&oflag))    oflag = 0;
    if (!parObj.getparint("cmp",&cmpflag))        cmpflag = 0;

    parObj.checkpars();

    /* absorb sector into scale */
    scale = scale/sector;

    /* get type and index values of output */
    type_out  = hdtype(key);
    index_out = getindex(key);
    if (oflag) {
        otype_out  = hdtype(okey);
        oindex_out = getindex(okey);    
    }
    if (cmpflag) {
        mxtype_out  = hdtype(mxkey);
        mxindex_out = getindex(mxkey);  
        mytype_out  = hdtype(mykey);
        myindex_out = getindex(mykey);
    }
    
    /* loop over traces */
    while (cs2su->getTrace(&tr)) {

        /* get header values */
        sx = (double) tr.sx;
        sy = (double) tr.sy;
        gx = (double) tr.gx;
        gy = (double) tr.gy;

        /* If tr.scalco not set, use 1 as the value */
        factor = (!tr.scalco) ? 1 : tr.scalco;

        /* factor < 0 means divide; factor > 0 means to multiply */
        if (factor < 0.0) factor = -1.0/factor;

	sx *= factor; sy *= factor;
	gx *= factor; gy *= factor;

        /* compute midpoints */
        mx = (double) (0.5*(sx + gx));
        my = (double) (0.5*(sy + gy));

        /* compute the azimuths */
        azimuth = computeAzimuth(sx,sy,mx,my,scale,az);

        /* set and put the new header value */
        setval(type_out,&val_out, azimuth);
        puthval(&tr, index_out, &val_out);

        /* set offset header word, if requested */
        if (oflag != 0) {
            int offsgn=0;

            /* calculate sign of offset */
            offsgn=SGN(MIN((sx-gx),(sy-gy)));
            
            /* calcluate offset */
            offset = sqrt((gx-sx)*(gx-sx) + (gy-sy)*(gy-sy));

            /* scale offset by scalco factor and set sign */
            offset *= offsgn;
            
            /* set and put the new header value */
            setval(otype_out,&val_out, offset);
            puthval(&tr, oindex_out, &val_out);
        }

        /* set cmp coord header words, if requested */
        if (cmpflag != 0) {
            cmpx = (sx + gx)/2/factor;
            cmpy = (sy + gy)/2/factor;
            
            /* set and put the new header values */
            setval(mxtype_out,&val_out, cmpx);
            puthval(&tr, mxindex_out, &val_out);

            setval(mytype_out,&val_out, cmpy);
            puthval(&tr, myindex_out, &val_out);
        }

        /* write traces */
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

/* Functions used internally */

#define DRADTODEG ((double) (180.0/PI))

double computeAzimuth(double sx, double sy, 
            double mx, double my, double scale, int az)
/************************************************************************
computeAzimuth - compute the azimuth of traces

*************************************************************************
Input:
sx      source x coordinate
sy      source y coordinate
mx      midpoint or geophone x coordinate
my      midpoint or geophone y coordinate
scale   scale factor
az      azimuth convention flag

Output:
        returns scaled azimuth in degrees

*************************************************************************
Notes:
Does computations as doubles. Outputs either as 0-179.999 or 
0-359.999 degrees, or in sectors of "sector" degrees in size.
Direction vector points from midpoint to source for az=1 and 
in opposite direction for az=-1. For az=0, reciprocity is assumed.
*************************************************************************
Author: CWP: John Stockwell, 28 Oct 1998
Modifications: Nils Maercklin put type conversion to function
               setval() and added az=-1 (az<0) option, Sep 2006
*************************************************************************/
{
    double dval_out=0.0;    /* output value (azimuth) */

    /* Compute the azimuth */
    if (az==0) { /* 0-179.999 convention */
        /* compute angle */
         dval_out = 90.0 - DRADTODEG*atan((sy-my)/(sx-mx)); 

        /* Make output go from 0-179.99999 degrees. */
        if (dval_out == 180.000) dval_out = 0.0;

    } 
    else if (az<0) { /* 0-359.999 from S to M */
        /* Compute angle */
        dval_out = DRADTODEG*atan2((mx-sx),(my-sy));

        /* Make output go from 0-359.99999 degrees. */
        if (dval_out < 0.) dval_out += 360.0;
    }
    else { /* 0-359.999 from M to S */
        /* Compute angle */
        dval_out = DRADTODEG*atan2((sx-mx),(sy-my));

        /* Make output go from 0-359.99999 degrees. */
        if (dval_out < 0.) dval_out += 360.0;
    }

    /* Return scaled azimuth */
    return dval_out*scale;
}


void setval( cwp_String type_out, Value *val_out, double dval_out)
/************************************************************************
setval - set value of a trace header word 

*************************************************************************
Input:
type_out    type of header word
dval_out    value of header word (double)

Output:
val_out     value according to type of header word

*************************************************************************
Credits: Einar Kajartansson, Jack K. Cohen, John Stockwell (sushw)
*************************************************************************/
{
    switch (*type_out) {
    case 's':
        throw cseis_geolib::csException("can't change char header word");
    break;
    case 'h':
        val_out->h = (short) NINT(dval_out);
    break;
    case 'u':
        val_out->u = (unsigned short) NINT(dval_out);
    break;
    case 'l':
        val_out->l = (long) NINT(dval_out);
    break;
    case 'v':
        val_out->v = (unsigned long) NINT(dval_out);
    break;
    case 'i':
        val_out->i = (int) NINT(dval_out);
    break;
    case 'p':
        val_out->p = (unsigned int) NINT(dval_out);
    break;
    case 'f':
        val_out->f = (float) dval_out;
    break;
    case 'd':
        val_out->d = (double) dval_out;
    break;
    default:
        throw cseis_geolib::csException("unknown type %s", type_out);
    break;
    }
}

} // END namespace
