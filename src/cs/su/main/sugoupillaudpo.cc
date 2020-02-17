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
std::string sdoc_sugoupillaudpo =
"									  "
" SUGOUPILLAUDPO - calculate Primaries-Only impulse response of a lossless"
"	      GOUPILLAUD medium for plane waves at normal incidence	"
"									  "
" sugoupillaudpo < stdin > stdout [optional parameters]		      "
"									    "
" Required parameters:							     "
"	none								"
"      									     "
" Optional parameters:						       "
"	l=1	   source layer number; 1 <= l <= tr.ns		  "
"		      Source is located at the top of layer l.		     "
"	k=1	   receiver layer number; 1 <= k			 "
"		      Receiver is located at the top of layer k.	    "
"	tmax	  number of output time-samples;			"
"		      default: long enough to capture all primaries	 "
"	pV=1	  flag for vector field seismogram		      "
"		      (displacement, velocity, acceleration);	       "
"		      =-1 for pressure seismogram.			  "
"	verbose=0     silent operation, =1 list warnings		    "
"									     "
" Input: Reflection coefficient series:				      "
"									    "
"			       impedance[i]-impedance[i+1]		   "
"		       r[i] = -----------------------------		  "
"			       impedance[i]+impedance[i+1]		   "
"									    "
"	r[0]= surface refl. coef. (as seen from above)		      "
"	r[n]= refl. coef. of the deepest interface			  "
"									     "
" Input file is to be in SU format, i.e., binary floats with a SU header.    "
"									    "
" Remarks:								   "
" 1. For vector fields, a buried source produces a spike of amplitude 1      "
" propagating downwards and a spike of amplitude -1 propagating upwards.     "
" A buried pressure source produces spikes of amplitude 1 both in the up-    "
" and downward directions.						   "
"    A surface source induces only a downgoing spike of amplitude 1 at the   "
" top of the first layer (both for vector and pressure fields).	      "
" 2. The sampling interval dt in the header of the input reflectivity file   "
" is interpreted as a two-way traveltime thicknes of the layers. The sampling"
" interval of the output seismogram is the same as that of the input file.   "
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sugoupillaudpo {


/* 
 * Credits:
 *	CWP: Albena Mateeva, April 2001.
 *
 */

/**************** end self doc ***********************************/

/* Functions used internally */

int imin( int c1, int c2 )
/* returns the smaller integer */
{
  if(c1<=c2)
    return(c1);
  else
    return(c2);
}

int imax( int c1, int c2 )
/* returns the larger integer */
{
  if(c1>=c2)
    return(c1);
  else
    return(c2);
}
/*********************************/

segy tr;

void* main_sugoupillaudpo( void* args )
{
  int n;	    /* number of subsurface interfaces		 */
  int l, k;	    /* source and receiver layers		      */
  int tmax;	    /* number of output samples			*/
  int pV;	    /* field-type flag				       */
  int verbose=0;    /* verbose flag				       */  
  float *r;	    /* input  reflectivity series		      */
  float *x;	    /* output seismogram			       */

  int i;	    /* loop counter      			       */
  int rmax;	 /* index of the deepest observable reflector       */
  int skl;	  /* sign of k-l				     */
  int n1, n2;       /* min and max of (k,l)			    */   
  float effsd=1;    /* effective source strength in downward direction */
  float transm1=1;  /* one-way transmission coefficient
		       from source to receiver			 */
  float transm2=1;  /* two-way transmission coefficient
		       between reflector and receiver		  */


  /* Initialize */
  cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
  cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
  cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
  int argc = suArgs->argc;
  char **argv = suArgs->argv;
  cseis_su::csSUGetPars parObj;

  void* retPtr = NULL;  /*  Dummy pointer for return statement  */
  su2cs->setSUDoc( sdoc_sugoupillaudpo );
  if( su2cs->isDocRequestOnly() ) return retPtr;
  parObj.initargs(argc, argv);

  try {  /* Try-catch block encompassing the main function body */

  
  /* Get subsurface model */
  if (!cs2su->getTrace(&tr)) throw cseis_geolib::csException("Can't get reflectivity");
  n = tr.ns - 1 ;

  /* Get source and receiver positions */
  if (!parObj.getparint("k",&k)) 	k=1;
  if (!parObj.getparint("l",&l))	l=1;

  /* compute some helpful quantities */
  if (k>l) skl=1;
  else if (k<l) skl=-1;
  else skl=0;
  n1=imin(k,l);
  n2=imax(k,l);

  /* get the rest of the parameters */
  if (!parObj.getparint("tmax",&tmax)) 
    tmax=(n2-n1+2+2*imax(imax(n1-1,n+1-n2),n+1-n1))/2;
  if (!parObj.getparint("pV",&pV))	pV=1;
  if (!parObj.getparint("verbose",&verbose))	verbose=0;
  
  /* Check parameters */
  if (n<0)  throw cseis_geolib::csException("The number of subsurface interfaces n cannot be negative.");
  if (n==0) warn("WARNING: model without subsurface reflectors!");
  if (k<1) throw cseis_geolib::csException("Receiver layer k must be >=1 (k=1 corresponds to surface seismogram).");
  if (l<1) throw cseis_geolib::csException("Source layer l must be >= 1 (l=1 corresponds to a surface source).");
  if (l>(n+1)) throw cseis_geolib::csException("The current version of the program requires l<=n+1.");
  if (tmax<0) throw cseis_geolib::csException("The number of the output time samples tmax cannot be negative.");
  if (! ( pV == 1 || pV == -1 ) ) throw cseis_geolib::csException("The field-type flag pV should be either 1 or -1.");
  
  /* Verbose */
  if (verbose) {
    warn("Source layer l=%d", l);
    warn("Number of layers n=%d", n);
    warn("Receiver layer k=%d", k);
    warn("Output time samples tmax=%d", tmax);
    warn("Field type flag pV=%d", pV);
  }
  
  /* Allocate Memory */
  r = ealloc1float(n+1);
  x = ealloc1float(2*tmax);
  
  /* Read and check reflectivity values */
  for(i=0; i<=n; ++i) {
    r[i] = tr.data[i];
    if(r[i]>1. || r[i]<-1.)
      throw cseis_geolib::csException("Invalid reflection coefficient encountered.");
  }
  
  /* Accounting for the field type (pressure/displacement) */
  for(i=0; i<=n; ++i)   r[i] *= pV;
  
  /* Initial zeroing of the seismogram */
  memset( (void *) x, 0, 2*tmax*FSIZE);

  /* Trivial case */
  if (n2-n1>=2*tmax)  throw cseis_geolib::csException("Seismogram length too small -- cannot observe any signal.");

  /* one-way transmission coef. between source and receiver */
  for(i=n1; i<n2 ; ++i) transm1 *= 1+skl*r[i];
  
  /* direct arrival */
  if (pV == -1)    x[n2-n1] = transm1;     /* r[i]=pV*r[i] already */
  else if (l == 1) x[n2-n1] = transm1;     /* downgoing 1 for pV=1 */
  else	     x[n2-n1] = skl*transm1; /* upgoing -1  for pV=1 */

  /* effective source in downward direction */
  if(!(l==1)) effsd=1+pV*r[l-1];
  
  /* index of the deepest observable reflector */
  i=(n2+n1-1)/2;
  rmax=imin(n,tmax-1+i);
  
  /* primary reflections from below the receiver */
  /* (and below the source if l>k) */
  for(i=n2; i<=rmax; ++i) 
    { 
      x[2*(i+1)-n2-n1] = effsd*transm1*r[i]*transm2; /* i<=n for sure */
      transm2 *= 1-r[i]*r[i]; 
    }
  
  /* primary reflections from above the receiver */
  /* (and above the source if l<k) */
  if (!(l==1)) 
	for(i=n1-1, transm2=1; i>=0 && n2-n1+2*(n1-1-i)<2*tmax; --i) 
      { 
	x[n2-n1+2*(n1-1-i)] += pV*transm1*r[i]*transm2; 
	transm2 *= 1-r[i]*r[i];
      } 

  /* output synthetic seismogram */
  i=(k-l)/2;
  if ( 2*i == k-l ){  
    for (i=0; i<tmax; ++i)  tr.data[i]=x[2*i];
    tr.delrt=0;
  }
  else{
    for (i=0; i<tmax; ++i)  tr.data[i]=x[2*i+1];
    tr.delrt=tr.dt/2000;
  }

  /* set new trace header */
  tr.ns=tmax;
  tr.trid=1;
	
  su2cs->putTrace(&tr);
  

  /* free allocated memmory */
  free1float(r);
  free1float(x);
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
