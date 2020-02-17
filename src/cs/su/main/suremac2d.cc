/* 
  Acoustic 2D Fourier method modeling with REM time integration
*/

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
extern "C" {
  #include <stdio.h>

}
extern "C" {
  #include <fcntl.h>

}
extern "C" {
  #include <math.h>

}
extern "C" {
  #include <stdlib.h>

}
extern "C" {
  #include <time.h>

}
#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
std::string sdoc_suremac2d =
"                                                                        "
" SUREMAC2D - Acoustic 2D Fourier method modeling with high accuracy     "
"             Rapid Expansion Method (REM) time integration              "
"                                                                        "
" suremac2d [parameters]                                                 "
"                                                                        "
" Required parameters:                                                   "
"                                                                        "
" opflag=     0: variable density wave equation                          "
"             1: constant density wave equation                          "
"             2: non-reflecting wave equation                            "
"                                                                        "
" nx=         number of grid points in horizontal direction              "
" nz=         number of grid points in vertical direction                "
" nt=         number of time samples                                     "
" dx=         spatial increment in horizontal direction                  "
" dz=         spatial increment in vertical direction                    "
" dt=         time sample interval in seconds                            "
" isx=        grid point # of horizontal source positions                "
" isz=        grid point # of vertical source positions                  "
"                                                                        "
" Optional parameters:                                                   "
" fx=0.0      first horizontal coordinate                                "
" fz=0.0      first vertical coordinate                                  "
" irx=        horizontal grid point # of vertical receiver lines         "
" irz=        vertical grid point # of horizontal receiver lines         "
" w=0.1       width of spatial source distribution (see notes)           "
" sflag=2     source time function                                       "
"             0: user supplied source function                           "
"             1: impulse (spike at t=0)                                  "
"             2: Ricker wavelet                                          "
" fmax=       maximum frequency of Ricker (default) wavelet              "
" amps=1.0    amplitudes of sources                                      "
" rfac=1.1    stability factor for expansion terms (see notes)           "
" prec=0      1: precompute Bessel coefficients b_k (see notes)          "
"             2: use precomputed Bessel coefficients b_k                 "
" fsflag=0    1: perform run with free surface b.c.                      "
" vmaxu=      user-defined maximum velocity                              "
" dtsnap=0.0  time interval in seconds of wave field snapshots           "
" iabso=1     apply absorbing boundary conditions (0: none)              "
" abso=0.1    damping parameter for absorbing boundaries                 "
" nbwx=20     horizontal width of absorbing boundary                     "
" nbwz=20     vertical width of absorbing boundary                       "
" verbose=0   1: show parameters used                                    "
"             2: print maximum amplitude at every expansion term         "
"                                                                        "
" velfile=vel          velocity filename                                 "
" densfile=dens        density filename                                  "
" sname=wavelet.su     user supplied source time function filename       "
" sepxname=sectx.su    x-direction pressure sections filename            "
" sepzname=sectz.su    z-direction pressure sections filename            "
" snpname=snap.su      pressure snapshot filename                        "
" jpfile=stderr        diagnostic output                                 "
"                                                                        "
" Notes:                                                                 "
"  0. The combination of the Fourier method with REM time integration    "
"     allows the computation of synthetic seismograms which are free     "
"     of numerical grid dispersion. REM has no restriction on the        "
"     time step size dt. The Fourier method requires at least two        "
"     grid points per shortest wavelength.                               "
"  1. nx and nz must be valid numbers for pfafft transform lengths.      "
"     nx and nz must be odd numbers (unless opflag=1). For valid         "
"     numbers see e.g. numbers in structure 'nctab' in source file       "
"     $CWPROOT/src/cwp/lib/pfafft.c.                                     "
"  2. Velocities (and densities) are stored as plain C style files       "
"     of floats where the fast dimension is along the z-direction.       "
"  3. Units must be consistent, e.g. m, s and m/s.                       "
"  4. A 20 grid points wide border at the sides and the bottom of        "
"     the modeling grid is used for sponge boundary conditions           "
"     (default: iabso=1).                                                "
"     Source and receiver lines should be placed some (e.g. 10) grid     "
"     points away from the absorbing boundaries in order to reduce       "
"     reflections due to obliquely incident wavefronts.                  "
"  5. Dominant frequency is about fmax/2 (sflag=2), absolute maximum     "
"     is delayed by 3/fmax from beginning of wavelet.                    "
"  6. If opflag!=1 the source should be not a spike in space; the        "
"     parameter w determines at which distance (in grid points) from     "
"     the source's center the Gaussian weight decays to 10 percent       "
"     of its maximum. w=2 may be a reasonable choice; however, the       "
"     waveform will be distorted.                                        "
"  7. Horizontal and vertical receiver line sections are written to      "
"     separate files. Each file can hold more than one line.             "
"  8. Parameter rfac may be enlarged if the modeling run becomes         "
"     unstable. This happens if the largest eigenvalue of the modeling   "
"     operator L is larger than estimated from the largest velocity.     "
"     In particular if using the variable density acoustic wave          "
"     equation the eigenvalues depend also on the density and it is      "
"     impossible to estimated the largest eigenvalue analytically.       "
"  9. Bessel coefficients can be precomputed (prec=1) and stored on      "
"     disk to save CPU time when several shots need to be run.           "
"     In this case computation of Bessel coefficients can be skipped     "
"     and read from disk file for reuse (prec=2).                        "
"     For reuse of Bessel coefficients the user may need to define       "
"     the overall maximum velocity (vmaxu).                              "
" 10. If snapshots are not required, a spike source (sflag=1) may be     "
"     applied and the resulting impulse response seismograms can be      "
"     convolved later with a desired wavelet.                            "
" 11. The free surface (fsflag=1) does not coincide with the first       "
"     vertical grid index (0). It appears to be half a grid spacing      "
"     above that position.                                               "
"                                                                        "
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suremac2d {

/**************** end self doc ***********************************/

/*
  Acoustic 2D Fourier method modeling with REM time integration

  Reference: 
  Kosloff, D., Fihlo A.Q., Tessmer, E. and Behle, A., 1989,
    Numerical solution of the acoustic and elastic wave equations by a
    new rapid expansion method, Geophysical Prospecting, 37, 383-394
  
 * Credits:
 *      University of Hamburg: Ekkehart Tessmer, October 2012
 */

FILE *snp, *sepx, *sepz, *vel, *dens, *jpfp, *sfp, *fpbes;

segy tri, tro, sno;

void* main_suremac2d( void* args )
{
  int i, k, l, il, it; /* loop indices */
  int indx, indz;
  int opflag;  /* selection of wave equation type */
  int sflag;   /* source function selection (Ricker/spike) */
  int fsflag;  /* flag for modeling with free surface b.c. */
  int prec;    /* flag for precomputing Bessel coefficients */
  int m;       /* number of Bessel coefficients used in expansion */
  int m0;
  int mm;      /* number of Bessel coefficients computed */
  int nx;      /* number of gridpoints in horizontal direction */
  int nz;      /* number of gridpoints in vertical direction */
  int nzpad;   /* number of gridpoints in vertical direction with zero-padding */
  int nt;      /* number of time samples */
  int nt0;
  int nbwx;    /* width of absorbing boundary in horizontal direction */
  int nbwz;    /* width of absorbing boundary in vertical direction */
  int nsnap;   /* number of snapshots */
  int nsnap0;
  int nsectx;  /* number of horizontal sections */
  int nsectz;  /* number of vertical sections */
  int nsourc;  /* number of sources */
  int iabso;   /* flag for using absorbing boundaries */
  int nwav=0;  /* number of samples for user source function */
  int verbose;
  int *irx=NULL; /* horizontal gridpoint indices of vertical receiver lines */
  int *irz=NULL; /* vertical gridpoint indices of horizontal receiver lines */
  int *isx;    /* horizontal gridpoint indices of sources */
  int *isz;    /* vertical gridpoint indices of sources */
  float abso;  /* damping value for absorbing boundaries */
  float dx;    /* grid spacing in horizontal direction */
  float dz;    /* grid spacing in vertical direction */
  float dt;    /* time step size or sample rate */
  float dt0;   
  float dtsnap; /* time increment between snapshots */
  float dtsnap0;
  float dtwave;/* time interval of user supplied wavelet (sflag=0) */
  float fx;    /* first horizontal coordinate */
  float fz;    /* first vertical coordinate */
  float t;     /* time */
  float tmax;  /* maximum time */
  float w;     /* width of source's spatial Gaussian distribution */
  float fmax;  /* maximum frequency of source time function */
  float vmax;  /* maximum subsurface velocity */
  float vmax0;
  float vmin;  /* minimum subsurface velocity */
  float vmaxu; /* user defined maximum velocity */
  float r;     /* largest eigenvalue of modeling operator L */
  float r0;
  float r2;    /* r^2 */
  float rfac;  /* factor for expansion terms */
  float amx;   /* maximum amplitude in numerical grid */
  float pi;
  float tmp;   /* temporary storage */
  float *amps;     /* amplitudes of sources */
  float *bwx;      /* weighting function values in absorbing zone (hor.) */
  float *bwz;      /* weighting function values in absorbing zone (vert.) */
  float *wave=NULL; /* user source function array */
  float **a1;       /* pressure wave field */
  float **a2;       /* pressure wave field */
  float **c2rho;    /* velocity structure c, c*c or c*c*rho (see opflag) */
  float **rhoinv=NULL; /* inverse of density */
  float **aux1;        /* auxiliary array */
  float **aux2=NULL;   /* auxiliary array */
  float **gbox=NULL;   /* Gaussian source box */
  float ***sectx=NULL; /* time sections along horizontal direction */
  float ***sectz=NULL; /* time sections along vertical direction */
  float ***snap=NULL;  /* snapshots */
  double **bessn, **bestr;
  char* jpfile;        /* file for information output */
  cwp_String velfile;  /* velocity filename */
  cwp_String densfile; /* density filename */
  cwp_String sname;    /* source function filename */
  cwp_String sepxname; /* x-direction pressure sections filename */
  cwp_String sepzname; /* z-direction pressure sections filename */
  cwp_String snpname;  /* pressure snapshot filename */

  /* function prototypes */
  void stru(int nx, int nz, float **c2rho, float **rhoinv, int opflag, float *vmin, float *vmax);
  void difx(float **vin, float **vout, float **bulk, int nx, int nz, float dx, int iadd);
  void difz(float **vin, float **vout, float **bulk, int nx, int nz, int nzpad, float dz, int iadd, int iload);
  void difx2(float **vin, float **vout, float **bulk, int nx, int nz, float dx, int iadd);
  void difz2(float **vin, float **vout, float **bulk, int nx, int nz, int nzpad, float dz, int iadd);
  void bessel_jn(double x, int n, double *bes);
  void intgr(float t, float r, float fmax, int mm, 
             float (*func) (float,float), double *sum);
  void intgru(float t, float r, float dt, int nt, int mm, 
             float *wave, double *sum);
  float fwave(float t, float fmax);
  void damp(int nbw, float abso, float *bw);
  void tbc2d(float *a, float *bwx, float *bwz, int nbwx, int nbwz,
	     int nx, int nz);
  void g2d(float w, int nx, int nz, float **g);

  /* defaults */
#define FX         0.0
#define FZ         0.0
#define IABSO      1
#define ABSO       0.1
#define NBWX      20
#define NBWZ      20
#define DTSNAP     0.0
#define SFLAG      2
#define NXB       21
#define NZB       21
#define W          0.1
#define RFAC       1.1
#define PREC       0
#define FSFLAG     0

  cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
  cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
  cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
  int argc = suArgs->argc;
  char **argv = suArgs->argv;
  cseis_su::csSUGetPars parObj;

  void* retPtr = NULL;  /*  Dummy pointer for return statement  */
  su2cs->setSUDoc( sdoc_suremac2d );
  if( su2cs->isDocRequestOnly() ) return retPtr;
  parObj.initargs(argc, argv);

  try {  /* Try-catch block encompassing the main function body */


  /* get parameters */
  if (!parObj.getparint ("opflag" , &opflag)) throw cseis_geolib::csException("opflag required!\n");
  if (!parObj.getparint ("nx" , &nx)) throw cseis_geolib::csException("nx required!\n");
  if (!parObj.getparint ("nz" , &nz)) throw cseis_geolib::csException("nz required!\n");
  if (!parObj.getparint ("nt" , &nt)) throw cseis_geolib::csException("nt required!\n");
  if (!parObj.getparfloat ("dx" , &dx)) throw cseis_geolib::csException("dx required!\n");
  if (!parObj.getparfloat ("dz" , &dz)) throw cseis_geolib::csException("dz required!\n");
  if (!parObj.getparfloat ("dt" , &dt)) throw cseis_geolib::csException("dt required!\n");
  
  nsourc = parObj.countparval("isx");
  if (nsourc != parObj.countparval("isz"))
    throw cseis_geolib::csException("isx and isz must have the same number of coordinates.\n");
  isx = ealloc1int(nsourc);
  isz = ealloc1int(nsourc);
  amps = ealloc1float(nsourc);
  for (l=0; l<nsourc; l++) amps[l] = 1.;
  if (!parObj.getparint ("isx" , isx)) throw cseis_geolib::csException("isx required!\n");
  if (!parObj.getparint ("isz" , isz)) throw cseis_geolib::csException("isz required!\n");
  
  /* get optional parameters */
  if (!parObj.getparfloat ("dtsnap" , &dtsnap)) dtsnap = DTSNAP;
  if (!parObj.getparint ("sflag" , &sflag)) sflag = SFLAG;
  if (!parObj.getparfloat ("fmax" , &fmax) && sflag==2) 
    throw cseis_geolib::csException("fmax required if Ricker wavelet selected!\n");
  if (!parObj.getparint ("prec" , &prec)) prec = PREC;
  if (!parObj.getparint ("fsflag" , &fsflag)) fsflag = FSFLAG;
  nsectx = parObj.countparval("irz");
  if (nsectx) {
    irz = ealloc1int(nsectx);
    parObj.getparint("irz", irz);
  }
  nsectz = parObj.countparval("irx");
  if (nsectz) {
    irx = ealloc1int(nsectz);
    parObj.getparint("irx", irx);
  }
  if (!parObj.getparint ("iabso" , &iabso)) iabso = IABSO;
  if (!parObj.getparfloat ("abso" , &abso)) abso = ABSO;
  if (!parObj.getparint ("verbose" , &verbose)) verbose=0;
  if (!parObj.getparint ("nbwx" , &nbwx)) nbwx = NBWX;
  if (!parObj.getparint ("nbwz" , &nbwz)) nbwz = NBWZ;
  if (!parObj.getparfloat ("w" , &w)) w = W;
  if (!parObj.getparfloat ("fx" , &fx)) fx = FX;
  if (!parObj.getparfloat ("fz" , &fz)) fz = FZ;
  if (!parObj.getparfloat ("vmaxu" , &vmaxu)) vmaxu = 0.;
  parObj.getparfloat("amps", amps);
  if (!parObj.getparfloat ("rfac" , &rfac)) rfac = RFAC;

  if (!parObj.getparstring("jpfile",&jpfile)) {
    jpfp = stderr;
  } else { 
    jpfp = fopen(jpfile,"w");
  }

  if (!parObj.getparstring("velfile",&velfile)) velfile = "vel"; 
  if (!parObj.getparstring("densfile",&densfile)) densfile = "dens"; 
  if (!parObj.getparstring("snpname",&snpname)) snpname = "snap.su"; 
  if (!parObj.getparstring("sepxname",&sepxname)) sepxname = "sectx.su"; 
  if (!parObj.getparstring("sepzname",&sepzname)) sepzname = "sectz.su"; 
  if (!parObj.getparstring("sname",&sname)) sname = "wavelet.su"; 

  nzpad = nz;
  if (fsflag == 1) {
    nzpad = npfa((int)(nz * 1.5));
    if (opflag == 1) {
      nzpad = npfa(nzpad);
    } else {
      while (!(nzpad % 2)) nzpad = npfa(++nzpad); /* force nzpad odd */
    }
  }
  
  if (verbose) {
    fprintf(jpfp,"opflag=%d\n",opflag);
    fprintf(jpfp,"nx=%d  dx=%f\n",nx,dx);
    fprintf(jpfp,"nz=%d  dz=%f\n",nz,dz);
    fprintf(jpfp,"fsflag=%d\n",fsflag);
    if (fsflag == 1) fprintf(jpfp,"nzpad=%d\n",nzpad);
    fprintf(jpfp,"nt=%d  dt=%f\n",nt,dt);
    fprintf(jpfp,"fx=%f  fz=%f\n",fx,fz);
    fprintf(jpfp,"w=%f\n",w);
    if (nsourc == 1) 
      fprintf(jpfp,"isx=%d  isz=%d\n",isx[0],isz[0]);
    else {
      fprintf(jpfp,"sources: nsourc=%d\n",nsourc);
      fprintf(jpfp,"  coordinates  amplitudes\n");
      for (il=0; il<nsourc; il++)
	fprintf(jpfp,"  (%4d,%4d)  %e\n",isx[il],isz[il],amps[il]);
    }
    fprintf(jpfp,"rfac=%e\n",rfac);
    fprintf(jpfp,"sflag=%d\n",sflag);
    fprintf(jpfp,"fmax=%e\n",fmax);
    fprintf(jpfp,"prec=%d\n",prec);

    if (nsectx) {
      fprintf(jpfp,"x-receiver lines: nsectx=%d\n",nsectx);
      fprintf(jpfp,"  z-coordinates:");
      for (il=0; il<nsectx; il++) fprintf(jpfp,"  %d",irz[il]);
      fprintf(jpfp,"\n");
    }
    if (nsectz) {
      fprintf(jpfp,"z-receiver lines: nsectx=%d\n",nsectx);
      fprintf(jpfp,"  x-coordinates:");
      for (il=0; il<nsectz; il++) fprintf(jpfp,"  %d",irx[il]);
      fprintf(jpfp,"\n");
    }
    if (vmaxu > 0.) fprintf(jpfp,"vmaxu=%f\n",vmaxu);
  }

  if (npfa(nx) != nx) throw cseis_geolib::csException("Error: Invalid nx !\n");
  if (! nx % 2 && opflag != 1) throw cseis_geolib::csException("Error: nx must be odd!\n");
  if (npfa(nz) != nz) throw cseis_geolib::csException("Error: Invalid nz !");
  if (! nz % 2 && opflag != 1) throw cseis_geolib::csException("Error: nz must be odd!\n");

  for (l=0; l<nsourc; l++) {
    if (isx[l] - (NXB+1)/2 < 0 && isx[l] + (NXB+1)/2 > nx-1) {
      throw cseis_geolib::csException("Error: (%d) source box beyond model's boundary in x-direction!\n",l);
    }
    if (isz[l] - (NZB+1)/2 < 0 && isz[l] + (NZB+1)/2 > nz-1) {
      throw cseis_geolib::csException("Error: (%d) source box beyond model's boundary in z-direction!\n",l);
    }
  }

  tmax = (nt-1) * dt;
  if (dtsnap < 0) throw cseis_geolib::csException("stopped: dtsnap must be >= 0!\n");
  if (dtsnap > 0.)
    nsnap = tmax / dtsnap + 1;
  else
    nsnap = 0;

  if (verbose) {
    if (nsectx) fprintf(jpfp,"sepxname=%s\n",sepxname);
    if (nsectz) fprintf(jpfp,"sepzname=%s\n",sepzname);
    fprintf(jpfp,"nsnap=%d\n",nsnap);
    if (nsnap) {
      fprintf(jpfp,"dtsnap=%e\n",dtsnap);
      fprintf(jpfp,"snpname=%s\n",snpname);
    }
  }
  
  if (prec == 2) 
    fprintf(jpfp,"sflag and fmax are ignored if precomputed Bessel coefficients are used (prec=2).\n");

  if (nt < 1) throw cseis_geolib::csException("stopped: nt must be > 0!\n");

  bwx    = ealloc1float(nbwx);
  bwz    = ealloc1float(nbwz);
  a1     = ealloc2float(nx,nz);
  a2     = ealloc2float(nx,nz);
  c2rho  = ealloc2float(nx,nz);
  aux1   = ealloc2float(nx,nzpad);
  gbox   = ealloc2float(NXB,NZB);

  if (opflag != 1) aux2 = ealloc2float(nx,nzpad);
  if (opflag == 0) rhoinv = ealloc2float(nx,nz);
  if (nsnap > 0) snap = ealloc3float(nx,nz,nsnap);
  if (nsectx > 0) sectx = ealloc3float(nx,nt,nsectx);
  if (nsectz > 0) sectz = ealloc3float(nz,nt,nsectz);

  dtwave = 0.;
  if (sflag == 0) {
    sfp = efopen(sname,"r");
    fgettr(sfp,&tri);
    nwav = tri.ns;
    dtwave = tri.dt * 1.e-6;
    wave = ealloc1float(nwav);
    for (it=0; it<nwav; it++) wave[it] = tri.data[it];
    efclose(sfp);
  }

  /* read velocity (and density) file */
  vel = efopen(velfile,"r");
  if (opflag == 0) dens = efopen(densfile,"r");

  stru(nx,nz,c2rho,rhoinv,opflag,&vmin,&vmax);

  /* change absorbing zone for no reflections at top */
  if (iabso) {
    for (k=0; k<nbwz; k++) {
      for (i=0; i<nx; i++) c2rho[nz-1-k][i] = c2rho[0][i];
      if (opflag == 0)
	for (i=0; i<nx; i++) rhoinv[nz-1-k][i] = rhoinv[0][i];
    }
  }

  efclose(vel);
  if (opflag == 0) efclose(dens);

  if (verbose) fprintf(jpfp,"max., min. velocities: %e, %e\n",vmax,vmin);
  
  /* setup aborbing boundaries */
  if (iabso) {
    damp(nbwx,abso,bwx);
    damp(nbwz,abso,bwz);
  }

  if (vmax == 0.) throw cseis_geolib::csException("stopped: vmax is zero!\n");
  
  if ((dx > dz ? dx : dz) > vmin/(2.*fmax))
    throw cseis_geolib::csException("stopped: fmax is too high for spatial sampling!\n");
  
  /* prepare Bessel coefficients convolved with wavelet */
  if (vmaxu > 0.) {
    if (vmax > vmaxu) throw cseis_geolib::csException("vmaxu is smaller than vmax!\n");
    vmax = vmaxu; 
  }
  pi = 4.*atan(1.);
  r = rfac * pi*sqrt(2.)*vmax/MIN(dx,dz); /* largest EV of operator L */
  r2 = r*r;

  m = tmax * r; /* highest index of terms in expansion */
  mm = m + 50;  /* added indices for better accuracy of high J_n(x) */

  if (verbose) fprintf(jpfp,"tmax*r=%e m=%d\n",tmax*r,m);
  bessn = ealloc2double(mm+1,nsnap);
  bestr = ealloc2double(mm+1,nt);

  /* Bessel-coefficients */
  switch(prec) {
  case 0:
  case 1:
    if (verbose) fprintf(jpfp,"computing Bessel coefficients ...");
    for (it=0; it<nsnap; it++) {
      t = it*dtsnap;
      if (sflag == 0)
	intgru(t,r,dtwave,nwav,mm,wave,bessn[it]);
      else if (sflag == 1) 
	bessel_jn(t*r,mm,bessn[it]);
      else 
	intgr(t,r,fmax,mm,fwave,bessn[it]);
    }
    
    for (it=0; it<nt; it++) {
      t = it*dt;
      if (sflag == 0)
	intgru(t,r,dtwave,nwav,mm,wave,bestr[it]);
      else if (sflag == 1) 
	bessel_jn(t*r,mm,bestr[it]);
      else 
	intgr(t,r,fmax,mm,fwave,bestr[it]);
    }
    if (verbose) fprintf(jpfp," done.\n");
    if (prec == 0) break;
    fpbes = efopen("b_k","w");
    efwrite(&vmax,sizeof(float),1,fpbes);
    efwrite(&r,sizeof(float),1,fpbes);
    efwrite(&m,sizeof(int),1,fpbes);
    efwrite(&nsnap,sizeof(int),1,fpbes);
    efwrite(&dtsnap,sizeof(float),1,fpbes);
    efwrite(&nt,sizeof(int),1,fpbes);
    efwrite(&dt,sizeof(float),1,fpbes);
    for (it=0; it<nsnap; it++) efwrite(bessn[it],sizeof(double),m,fpbes);
    for (it=0; it<nt; it++) efwrite(bestr[it],sizeof(double),m,fpbes);
    fclose(fpbes);
    break;
  case 2:
    if (!(fpbes = fopen("b_k","r")))
	throw cseis_geolib::csException("Bessel coefficients must be precomputed (prec=1).\n");
    efread(&vmax0,sizeof(float),1,fpbes);
    efread(&r0,sizeof(float),1,fpbes);
    efread(&m0,sizeof(int),1,fpbes);
    efread(&nsnap0,sizeof(int),1,fpbes);
    efread(&dtsnap0,sizeof(float),1,fpbes);
    efread(&nt0,sizeof(int),1,fpbes);
    efread(&dt0,sizeof(float),1,fpbes);
    if (vmax != vmax0 || r != r0 || m != m0 || nsnap != nsnap0 ||
	dtsnap != dtsnap0 || nt != nt0 || dt != dt0 )
      throw cseis_geolib::csException("Wrong set of Bessel coefficients in file b_k!\n");
    for (it=0; it<nsnap; it++) efread(bessn[it],sizeof(double),m,fpbes);
    for (it=0; it<nt; it++) efread(bestr[it],sizeof(double),m,fpbes);
    fclose(fpbes);
  }

  if (verbose) fprintf(jpfp,"starting modeling.\n");

  for (k=0; k<nz; k++) {
    for (i=0; i<nx; i++) a1[k][i] = a2[k][i] = 0.;
  }

  /* Gaussian source box */
  g2d(w, NXB, NZB, gbox);

  /* source distribution */
  for (l=0; l<nsourc; l++) {
    for (k=0; k<NZB; k++) {
      indz = isz[l]+k-NZB/2;
      for (i=0; i<NXB; i++) {
	indx = isx[l]+i-NXB/2;
	if (gbox[k][i] > 1.e-2) {
	  if (indx<0 || indx>nx-1 || indz<0 || indz>nz-1)
	    throw cseis_geolib::csException("source box beyond grid boundary!\n");
	  a1[indz][indx] = gbox[k][i] * amps[l];
	}
      }
    }
  }

  /* initialize first two expansion terms Q_1 and Q_3 */
  switch(opflag) {
  case 0:
    difx(a1,aux1,rhoinv,nx,nz,dx,2);
    difx(aux1,aux1,NULL,nx,nz,dx,0);

    difz(a1,aux2,rhoinv,nx,nz,nzpad,dz,2,0);
    difz(aux2,aux1,NULL,nx,nz,nzpad,dz,1,1);
    
    for (k=0; k<nz; k++) {
      for (i=0; i<nx; i++) {
	aux1[k][i] = c2rho[k][i]*aux1[k][i]/r2;
	a2[k][i] = 4.*aux1[k][i] + 3.*a1[k][i];
      }
    }
    break;
  case 1:
    difx2(a1,aux1,NULL,nx,nz,dx,0);
    difz2(a1,aux1,NULL,nx,nz,nzpad,dz,1);
    
    for (k=0; k<nz; k++) {
      for (i=0; i<nx; i++) {
	aux1[k][i] = c2rho[k][i]*aux1[k][i]/r2;
	a2[k][i] = 4.*aux1[k][i] + 3.*a1[k][i];
      }
    }
    break;
  case 2:
    difx(a1,aux1,c2rho,nx,nz,dx,2);
    difx(aux1,aux1,NULL,nx,nz,dx,0);

    difz(a1,aux2,c2rho,nx,nz,nzpad,dz,2,0);
    difz(aux2,aux1,NULL,nx,nz,nzpad,dz,1,1);
    
    for (k=0; k<nz; k++) {
      for (i=0; i<nx; i++) {
	aux1[k][i] = c2rho[k][i]*aux1[k][i]/r2;
	a2[k][i] = 4.*aux1[k][i] + 3.*a1[k][i];
      }
    }
    break;
  default:
    throw cseis_geolib::csException("no such opflag!\n");
    break;
  } /* end switch */

  if (nsnap > 0) {
    for (it=0; it<nsnap; it++) { /* first two terms for snapshots */
      for (k=0; k<nz; k++) {
	for (i=0; i<nx; i++) {
	  snap[it][k][i] = 2.*a1[k][i]*bessn[it][1] + 2.*a2[k][i]*bessn[it][3];
	}
      }
    }
  }
  
  if (nsectx > 0) {
    for (il=0; il<nsectx; il++) {
      for (it=0; it<nt; it++) { /* first two terms for x-sections */
	for (i=0; i<nx; i++) {
	  sectx[il][it][i] = 2.*a1[irz[il]][i]*bestr[it][1] + 2.*a2[irz[il]][i]*bestr[it][3]; 
	}
      }
    }
  }
  
  if (nsectz > 0) {
    for (il=0; il<nsectz; il++) {
      for (it=0; it<nt; it++) { /* first two terms for z-sections */
	for (k=0; k<nz; k++) {
	  sectz[il][it][k] = 2.*a1[k][irx[il]]*bestr[it][1] + 2.*a2[k][irx[il]]*bestr[it][3]; 
	}
      }
    }
  }

  if (iabso) {
    tbc2d(a1[0],bwx,bwz,nbwx,nbwz,nx,nz);
    tbc2d(a2[0],bwx,bwz,nbwx,nbwz,nx,nz);
  }
  
  /* calculate all other odd expansion terms Q_2n+1 */
  for (l=2; l<=m/2; l++) {
    if (verbose == 2) {
      amx = 0.;
      for (k=0; k<nz; k++) {
	for (i=0; i<nx; i++) {
	  tmp = fabs(a1[k][i]);
	  amx = (amx > tmp ? amx : tmp);
	}
      }
      fprintf(jpfp,"l=%4d  max= %e\n",l,amx);
    }

    switch(opflag) {
    case 0:
      difx(a2,aux1,rhoinv,nx,nz,dx,2);
      difx(aux1,aux1,NULL,nx,nz,dx,0);

      difz(a2,aux2,rhoinv,nx,nz,nzpad,dz,2,0);
      difz(aux2,aux1,NULL,nx,nz,nzpad,dz,1,1);

      for (k=0; k<nz; k++) {
	for (i=0; i<nx; i++) {
	  aux1[k][i] = c2rho[k][i]*aux1[k][i]/r2;
	  aux1[k][i] = 4.*aux1[k][i] + 2.*a2[k][i] - a1[k][i];
	  
	  a1[k][i] = a2[k][i];
	  a2[k][i] = aux1[k][i];
	}
      }
      break;
    case 1:
      difx2(a2,aux1,NULL,nx,nz,dx,0);
      difz2(a2,aux1,NULL,nx,nz,nzpad,dz,1);

      for (k=0; k<nz; k++) {
	for (i=0; i<nx; i++) {
	  aux1[k][i] = c2rho[k][i]*aux1[k][i]/r2;
	  aux1[k][i] = 4.*aux1[k][i] + 2.*a2[k][i] - a1[k][i];
	  
	  a1[k][i] = a2[k][i];
	  a2[k][i] = aux1[k][i];
	}
      }
      break;
    case 2:
      difx(a2,aux1,c2rho,nx,nz,dx,2);
      difx(aux1,aux1,NULL,nx,nz,dx,0);

      difz(a2,aux2,c2rho,nx,nz,nzpad,dz,2,0);
      difz(aux2,aux1,NULL,nx,nz,nzpad,dz,1,1);

      for (k=0; k<nz; k++) {
	for (i=0; i<nx; i++) {
	  aux1[k][i] = c2rho[k][i]*aux1[k][i]/r2;
	  aux1[k][i] = 4.*aux1[k][i] + 2.*a2[k][i] - a1[k][i];
	  
	  a1[k][i] = a2[k][i];
	  a2[k][i] = aux1[k][i];
	}
      }
      break;
    default:
      throw cseis_geolib::csException("no such opflag!\n");
      break;
    } /* end switch */

    if (nsnap > 0) {
      for (it=0; it<nsnap; it++) {
	for (k=0; k<nz; k++) { /* add term for snapshots */
	  for (i=0; i<nx; i++) {
	    snap[it][k][i] += 2.*a2[k][i]*bessn[it][2*l+1];
	  }
	}
      }
    }
      
    if (nsectx > 0) {	/* add term for x-sections */
      for (il=0; il<nsectx; il++) {
	for (it=0; it<nt; it++) {
	  for (i=0; i<nx; i++) {
	    sectx[il][it][i] += 2.*a2[irz[il]][i]*bestr[it][2*l+1];
	  }
	}
      }
    }
      
    if (nsectz > 0) { /* add term for z-sections */
      for (il=0; il<nsectz; il++) {
	for (it=0; it<nt; it++) {
	  for (k=0; k<nz; k++) {
	    sectz[il][it][k] += 2.*a2[k][irx[il]]*bestr[it][2*l+1];
	  }
	}
      }
    }

    if (iabso) {
      tbc2d(a1[0],bwx,bwz,nbwx,nbwz,nx,nz);
      tbc2d(a2[0],bwx,bwz,nbwx,nbwz,nx,nz);
    }

  } /* end loop over expansion terms */
    
  /* --- output results to files --- */
  /* common trace headers */
  tro.trid = 1;
  tro.scalco = -10;
  tro.sx = (int)((fx+isx[0]*dx)*10);
  tro.ns = nt;
  tro.dt = (int)(dt*1000000+0.5);

  /* write x-direction sections */
  if (nsectx) {
    sepx = efopen(sepxname,"w");
    for (il=0; il<nsectx; il++) {
      tro.fldr = il+1;
      for (i=0; i<nx; i++) {
	tro.tracl = i;
	tro.tracf = i;
	tro.gx = (int)((fx+i*dx)*10);
	for (it=0; it<nt; it++) tro.data[it] = sectx[il][it][i];
	fvsu2cs->putTrace(sepx,&tro);
      }
    }
    efclose(sepx);
  }

  /* write z-direction sections */
  if (nsectz) {
    sepz = efopen(sepzname,"w");
    for (il=0; il<nsectz; il++) {
      tro.fldr = il+1;
      for (k=0; k<nz; k++) {
	tro.tracl = k;
	tro.tracf = k;
	tro.gx = (int)((fz+k*dz)*10);
	for (it=0; it<nt; it++) tro.data[it] = sectz[il][it][k];
	fvsu2cs->putTrace(sepz,&tro);
      }
    }
    efclose(sepz);
  }

  /* write snapshots */
  if (nsnap) {
    snp = efopen(snpname,"w");
    sno.trid = 130; /* depth range */
    sno.scalco = -10;
    sno.sx = (int)((fx+isx[0]*dx)*10);
    sno.ns = nz;
    sno.dt = dz*10;
    sno.d1 = dz;
    sno.d2 = dx;
    for (l=0; l<nsnap; l++) {
      sno.fldr = l+1;
      for (i=0; i<nx; i++) {
	sno.tracl = i;
	sno.tracf = i;
	sno.gx = (int)((fx+i*dx)*10);
	sno.offset = (sno.gx-sno.sx)/10.;
	for (k=0; k<nz; k++) sno.data[k] = snap[l][k][i];
	fvsu2cs->putTrace(snp,&sno);
      }
    }
    efclose(snp);
  }

  efclose(jpfp);
  if (verbose) fprintf(jpfp,"finished modeling.");

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
/**********************************************************************/
void stru(int nx, int nz, float **c2rho, float **rhoinv, int opflag, float *vmin, float *vmax)
/*
  read structure velocities (and densities) file

  input:
    nx:      number of gridpoints in horizontal direction
    nz:      number of gridpoints in vertical direction
    opflag:  flag for wave equation type 

  output:
    c2rho:   array of c*c*rho or c*c or c values
    rhoinv:  array of 1/rho
    vmin:    minimum velocity of grid
    vmax:    maximum velocity of grid
*/
  {
  int i, k;
  float tmp;
  float *bufz;

  bufz = ealloc1float(nz);

  for (i=0; i<nx; i++) {
    efread(bufz,sizeof(float),nz,vel); 
    for (k=0; k<nz; k++) c2rho[k][i] = bufz[k];
  }

  switch(opflag) {
  case 0:
    /* variable density case */
    for (i=0; i<nx; i++) {
      efread(bufz,sizeof(float),nz,dens); 
      for (k=0; k<nz; k++) {
	rhoinv[k][i] = 1./bufz[k];
	c2rho[k][i] = c2rho[k][i]*c2rho[k][i]*bufz[k];
      }
    }
    break;
  case 1:
    /* constant density case */
    for (k=0; k<nz; k++) {
      for (i=0; i<nx; i++) {
	c2rho[k][i] = c2rho[k][i]*c2rho[k][i];
      }
    }
    break;
  case 2:
    /* non-reflecting wave equation */
    break;
  }

  /* find maximum and minimum velocity */
  *vmax = 0.;
  *vmin = 1.e20;
  switch(opflag) {
  case 0:
    for (k=0; k<nz; k++) {
      for (i=0; i<nx; i++) {
	tmp = sqrt(c2rho[k][i]*rhoinv[k][i]);
	*vmax = (tmp > *vmax ? tmp : *vmax);
	*vmin = (tmp < *vmin ? tmp : *vmin);
      }
    }
    break;
  case 1:
    for (k=0; k<nz; k++) {
      for (i=0; i<nx; i++) {
	tmp = sqrt(c2rho[k][i]);
	*vmax = (tmp > *vmax ? tmp : *vmax);
	*vmin = (tmp < *vmin ? tmp : *vmin);
      }
    }
    break;
  case 2:
    for (k=0; k<nz; k++) {
      for (i=0; i<nx; i++) {
	tmp = c2rho[k][i];
	*vmax = (tmp > *vmax ? tmp : *vmax);
	*vmin = (tmp < *vmin ? tmp : *vmin);
      }
    }
    break;
  }
  free1float(bufz);
}
/**********************************************************************/
void g2d(float w, int nx, int nz, float **g)
/*
  2-D Gaussian weights for spatial distribution of source 

  input:
    w:     width in grid points where Gaussian is one tenth of its maximum
    nx:    number of gridpoints in horizontal direction
    nz:    number of gridpoints in vertical direction

  output:
    g:     weights for spatial distribution of source
*/
{
  int i, i0, k, k0;
  float pi, sigma, fac, sum, r2, arg;

  pi = 4.*atan(1.);

  i0 = nx/2+1;
  k0 = nz/2+1;

  if (w < 0.1) w = 0.1; /* spatial spike */

  /* w is width where Gaussian is one tenth of its maximum */
  sigma = sqrt(-w*w/(2.*log(0.1)));

  fac = 1./(sigma*sqrt(2.*pi));

  sum = 0.;
  for (k=0; k<nz; k++) {
    for (i=0; i<nx; i++) {
      r2 = (i-i0)*(i-i0)+(k-k0)*(k-k0);
      arg = r2/(2.*sigma*sigma);
      g[k][i] = 0.;
      if (arg < 80.) g[k][i] = fac*exp(-arg);
      sum += g[k][i];
    }
  }
  /* normalization */
  for (k=0; k<nz; k++) {
    for (i=0; i<nx; i++) g[k][i] /= sum;
  }
}
/**********************************************************************/
void difx(float **vin, float **vout, float **bulk, int nx, int nz, float dx, 
	  int iadd)
/* 
   compute 1st derivative along x-direction 

   input:
     vin:    2D input array
     bulk:   2D array for multiplication with output array
     nx:     length of array in x-direction
     nz:     length of array in z-direction
     dx:     grid spacing in x-direction
     iadd:   flag

   output:
     vout:   2D output array

     iadd=0: set result in vout
     iadd=1: add result to vout
     iadd=2: set result multiplied by bulk in vout
     iadd=3: add result multiplied by bulk to vout
*/
   
{
  static int i, k;
  float  tmp;
  static float *a, *rkx;
  static int icall=0, ieo=0, m1=-1, p1=1 ;
  static int num;

  void rk(float *ak, int n, float d, int ind);
 
  if (icall == 0 ) {
    icall = 1;
    ieo = nz % 2;
    num = ((nz + 1) / 2) * 2;

    a   = ealloc1float(2*nx);
    rkx = ealloc1float(nx);
    rk(rkx,nx,dx,1);
  }

  for (k=0; k<num; k+=2) {
	
    /* Load data from vin */
    for (i=0; i<nx; i++) {
      a[2*i] = vin[k][i];
      if (!ieo || !(k==num-2)) {
	a[2*i+1] = vin[k+1][i];
      }
      else {
	a[2*i+1] = 0.;
      }
    }

    /* Perform derivative */
    pfacc(m1, nx, (complex*)a);
	
    for (i=0; i<nx; i++) {
      tmp      =  a[2*i]   * rkx[i];
      a[2*i]   = -a[2*i+1] * rkx[i];
      a[2*i+1] =  tmp;
    }
	
    pfacc(p1, nx, (complex*)a);

    if (iadd == 0) {

      /* Set data into vout */
      for (i=0; i<nx; i++) {
	vout[k][i] = a[2*i];
	if (!ieo || !(k==num-2)) vout[k+1][i] = a[2*i+1];
      }
    } else if (iadd == 1) {

      /* Add data into vout */
      for (i=0; i<nx; i++) {
	vout[k][i] += a[2*i];
	if (!ieo || !(k==num-2)) vout[k+1][i] += a[2*i+1];
      }
    } else if (iadd == 2) {

      /* Set data multiplied by bulk into vout */
      for (i=0; i<nx; i++) {
	vout[k][i] = a[2*i] * bulk[k][i];
	if (!ieo || !(k==num-2)) vout[k+1][i] = a[2*i+1] * bulk[k+1][i];
      }
    } else if (iadd == 3) {

      /* Add data multiplied by bulk into vout */
      for (i=0; i<nx; i++) {
	vout[k][i] += a[2*i] * bulk[k][i];;
	if (!ieo || !(k==num-2)) vout[k+1][i] += a[2*i+1] * bulk[k+1][i];
      }
    }
  }
}
/**********************************************************************/
void difx2(float **vin, float **vout, float **bulk, int nx, int nz, float dx, 
   int iadd)
/* 
   compute 2nd derivative along x-direction 

   input:
     vin:    2D input array
     bulk:   2D array for multiplication with output array
     nx:     length of array in x-direction
     nz:     length of array in z-direction
     dx:     grid spacing in x-direction
     iadd:   flag

   output:
     vout:   2D output array

   iadd=0: set result in vout
   iadd=1: add result to vout
   iadd=2: set result multiplied by bulk in vout
   iadd=3: add result multiplied by bulk to vout
*/

{
  static int i, k;
  static float *a, *rkx;
  static int icall=0, ieo=0, m1=-1, p1=1 ;
  static int num;

  void rk(float *ak, int n, float d, int ind);
 
  if (icall == 0 ) {
    icall = 1;
    ieo = nz % 2;
    num = ((nz + 1) / 2) * 2;

    a   = ealloc1float(2*nx);
    rkx = ealloc1float(nx);
    rk(rkx,nx,dx,2);
  }

  for (k=0; k<num; k+=2) {
	
    /* Load data from vin */
    for (i=0; i<nx; i++) {
      a[2*i] = vin[k][i];
      if (!ieo || !(k==num-2)) {
	a[2*i+1] = vin[k+1][i];
      }
      else {
	a[2*i+1] = 0.;
      }
    }

    /* Perform derivative */
    pfacc(m1, nx, (complex*)a);
	
    for (i=0; i<nx; i++) {
      a[2*i]   *= rkx[i];
      a[2*i+1] *= rkx[i];
    }
	
    pfacc(p1, nx, (complex*)a);

    if (iadd == 0) {

      /* Set data into vout */
      for (i=0; i<nx; i++) {
	vout[k][i] = a[2*i];
	if (!ieo || !(k==num-2)) vout[k+1][i] = a[2*i+1];
      }
    } else if (iadd == 1) {

      /* Add data into vout */
      for (i=0; i<nx; i++) {
	vout[k][i] += a[2*i];
	if (!ieo || !(k==num-2)) vout[k+1][i] += a[2*i+1];
      }
    } else if (iadd == 2) {

      /* Set data multiplied by bulk into vout */
      for (i=0; i<nx; i++) {
	vout[k][i] = a[2*i] * bulk[k][i];
	if (!ieo || !(k==num-2)) vout[k+1][i] = a[2*i+1] * bulk[k+1][i];
      }
    } else if (iadd == 3) {

      /* Add data multiplied by bulk into vout */
      for (i=0; i<nx; i++) {
	vout[k][i] += a[2*i] * bulk[k][i];;
	if (!ieo || !(k==num-2)) vout[k+1][i] += a[2*i+1] * bulk[k+1][i];
      }
    }
  }
}
/**********************************************************************/
void difz(float **vin, float **vout, float **bulk, int nx, int nz, int nzpad, float dz,
	  int iadd, int iload)
/* 
   compute 1st derivative along z-direction 

   input:
     vin:    2D input array
     bulk:   2D array for multiplication with output array
     nx:     length of array in x-direction
     nz:     length of array in z-direction
     dz:     grid spacing in z-direction
     iadd:   flag

   output:
     vout:   2D output array

   iadd=0: set result in vout
   iadd=1: add result to vout
   iadd=2: set result multiplied by bulk in vout
   iadd=3: add result multiplied by bulk to vout

   iload=0: load nz rows and zero-pad, unload nzpad rows
   iload=1: load nzpad rows, unload nz rows
*/

{
  static int i, k, nzunl;
  float  tmp;
  static float *a, *rkz;
  static int icall=0, ieo=0, m1=-1, p1=1 ;
  static int num;

  void rk(float *ak, int n, float d, int ind);

  if (icall == 0 ) {
    icall = 1;
    ieo = nx % 2;
    num = ((nx + 1) / 2) * 2;

    a   = ealloc1float(2*nzpad);
    rkz = ealloc1float(nzpad);
    rk(rkz,nzpad,dz,1);
  }

  for (i=0; i<num; i+=2) {

    if (iload == 0) {
      /* Load data from vin */
      for (k=0; k<nz; k++) {
	a[2*k] = vin[k][i];
	if (!ieo || !(i==num-2)) {
	  a[2*k+1] = vin[k][i+1];
	}
	else {
	  a[2*k+1] = 0.;
	}
      }
      /* zero padding for free surface */
      for (k=nz; k<nzpad; k++) {
	a[2*k]   = 0.;
	a[2*k+1] = 0.;
      }      
      nzunl = nzpad;
    } else {
      for (k=0; k<nzpad; k++) {
	a[2*k] = vin[k][i];
	if (!ieo || !(i==num-2)) {
	  a[2*k+1] = vin[k][i+1];
	}
	else {
	  a[2*k+1] = 0.;
	}
      }
      nzunl = nz;
    }

    /* Perform derivative */
    pfacc(m1, nzpad, (complex*)a);
	
    for (k=0; k<nzpad; k++) {
      tmp      =  a[2*k]   * rkz[k];
      a[2*k]   = -a[2*k+1] * rkz[k];
      a[2*k+1] =  tmp;
    }
	
    pfacc(p1, nzpad, (complex*)a);

    if (iadd == 0) {
	
      /* Set data into vout */
      for (k=0; k<nzunl; k++) {
	vout[k][i] = a[2*k];
	if (!ieo || !(i==num-2)) vout[k][i+1] = a[2*k+1];
      }
    } else if (iadd == 1) {

      /* Add data into vout */
      for (k=0; k<nzunl; k++) {
	vout[k][i] += a[2*k];
	if (!ieo || !(i==num-2)) vout[k][i+1] += a[2*k+1];
      }
    } else if (iadd == 2) {

      /* Set data multiplied by bulk into vout */
      for (k=0; k<nz; k++) {
	vout[k][i] = a[2*k] * bulk[k][i];
	if (!ieo || !(i==num-2)) vout[k][i+1] = a[2*k+1] * bulk[k][i+1];
      }
      if (nzunl == nzpad) {
	for (k=nz; k<nzpad; k++) {
	  vout[k][i] = a[2*k] * bulk[0][i];
	  if (!ieo || !(i==num-2)) vout[k][i+1] = a[2*k+1] * bulk[0][i+1];
	}
      }
    } else if (iadd == 3) {

      /* Add data multiplied by bulk into vout */
      for (k=0; k<nz; k++) {
	vout[k][i] += a[2*k] * bulk[k][i];
	if (!ieo || !(i==num-2)) vout[k][i+1] += a[2*k+1] * bulk[k][i+1];
      }
      if (nzunl == nzpad) {
	for (k=nz; k<nzpad; k++) {
	  vout[k][i] += a[2*k] * bulk[0][i];
	  if (!ieo || !(i==num-2)) vout[k][i+1] += a[2*k+1] * bulk[0][i+1];
	}
      }
    }
  }
}
/**********************************************************************/
void difz2(float **vin,float **vout, float **bulk, int nx, int nz, int nzpad, float dz,
	   int iadd)
/* 
   compute 2ndq derivative along z-direction 

   input:
     vin:    2D input array
     bulk:   2D array for multiplication with output array
     nx:     length of array in x-direction
     nz:     length of array in z-direction
     dz:     grid spacing in z-direction
     iadd:   flag

   output:
     vout:   2D output array

   iadd=0: set result in vout
   iadd=1: add result to vout
   iadd=2: set result multiplied by bulk in vout
   iadd=3: add result multiplied by bulk to vout
*/

{
  static int i, k;
  static float *a, *rkz;
  static int icall=0, ieo=0, m1=-1, p1=1 ;
  static int num;

  void rk(float *ak, int n, float d, int ind);

  if (icall == 0 ) {
    icall = 1;
    ieo = nx % 2;
    num = ((nx + 1) / 2) * 2;

    a   = ealloc1float(2*nzpad);
    rkz = ealloc1float(nzpad);
    rk(rkz,nzpad,dz,2);
  }

  for (i=0; i<num; i+=2) {
	
    /* Load data from vin */
    for (k=0; k<nz; k++) {
      a[2*k] = vin[k][i];
      if (!ieo || !(i==num-2)) {
	a[2*k+1] = vin[k][i+1];
      }
      else {
	a[2*k+1] = 0.;
      }
    }

    /* zero padding for free surface */
    for (k=nz; k<nzpad; k++) {
      a[2*k]   = 0.;
      a[2*k+1] = 0.;
    }      

    /* Perform derivative */
    pfacc(m1, nzpad, (complex*)a);
	
    for (k=0; k<nzpad; k++) {
      a[2*k]   *= rkz[k];
      a[2*k+1] *= rkz[k];
    }
	
    pfacc(p1, nzpad, (complex*)a);

    if (iadd == 0) {
	
      /* Set data into vout */
      for (k=0; k<nzpad; k++) {
	vout[k][i] = a[2*k];
	if (!ieo || !(i==num-2)) vout[k][i+1] = a[2*k+1];
      }
    } else if (iadd == 1) {

      /* Add data into vout */
      for (k=0; k<nzpad; k++) {
	vout[k][i] += a[2*k];
	if (!ieo || !(i==num-2)) vout[k][i+1] += a[2*k+1];
      }
    } else if (iadd == 2) {

      /* Set data multiplied by bulk into vout */
      for (k=0; k<nz; k++) {
	vout[k][i] = a[2*k] * bulk[k][i];
	if (!ieo || !(i==num-2)) vout[k][i+1] = a[2*k+1] * bulk[k][i+1];
      }
      for (k=nz; k<nzpad; k++) {
	vout[k][i] = a[2*k] * bulk[0][i];
	if (!ieo || !(i==num-2)) vout[k][i+1] = a[2*k+1] * bulk[0][i+1];
      }
    } else if (iadd == 3) {

      /* Add data multiplied by bulk into vout */
      for (k=0; k<nz; k++) {
	vout[k][i] += a[2*k] * bulk[k][i];
	if (!ieo || !(i==num-2)) vout[k][i+1] += a[2*k+1] * bulk[k][i+1];
      }
      for (k=nz; k<nzpad; k++) {
	vout[k][i] += a[2*k] * bulk[0][i];
	if (!ieo || !(i==num-2)) vout[k][i+1] += a[2*k+1] * bulk[0][i+1];
      }
    }
  }
}
 
/**********************************************************************/
void rk(float *ak, int n, float d, int ind)
/*
  wave-numbers for Fourier transform

  input:
    n:      length of Fourier transform
    d:      grid spacing
    ind:    order of derivative

  output:
    ak:     array of wave numbers
*/
{
  int i,isign,n2;
  float pi, dn, c; 
  
  pi = 4.*atan(1.);
  isign = (int)(-pow((-1.),(float)ind));
  n2 = n/2;
  dn = 2.*pi/(n*d);
  
  for (i=0; i<n; i++) {
    c=dn*i;
    if (i > n2) c = -dn*(n-i);
    ak[i] = pow(c,(float)ind)*isign/n;
  }
}
/**********************************************************************/
void bessel_jn(double x, int n, double *bes)
/*
  compute Bessel function values J0(x), J1(x), ...

  input:
    x:    argument
    n:    highest order of Bessel function

  output:
    bes:  function values of n orders of Bessel functions at argument x

  Notes:  1.  n must be larger than x
          2.  for accuracy of higher order Bessel function  choose n 
              larger (e.g. +50) than needed
*/
{
  int i, j;
  double sum;
  double eps=1.0e-130;

  if (x > (double)n)
    throw cseis_geolib::csException("bessel_jn: n(%d) must be larger than x(%e)!\n",n,x);

  if (fabs(x) > eps) {

    bes[n-1] = 0.;
    bes[n-2] = 1.;
    
    sum = 0.;
    for (i=n-3; i>=0; i--) {
      bes[i] = (2*(i+1)/x)*bes[i+1] - bes[i+2];
      
      if (i%2 == 0) {
	sum = sum + 2. * bes[i];
	if (fabs(sum) > 1.e30) {
	  for (j=i; j<n; j++) bes[j] /= sum;
	  sum = 1.;
	}
      }
    }
    
    sum = bes[0];
    for (i=2; i<n; i+=2) sum += 2. * bes[i];
    for (i=0; i<n; i++) bes[i] /= sum;
  } else { /* if fabs(x) very small assume x=0. */
    bes[0] = 1.;
    for (i=1; i<n; i++) bes[i] = 0.;
  }
}

/**********************************************************************/
void intgr(float t, float r, float fmax, int mm, float (*func)(float,float), double* sum)
/*
  integrate J_n(tR) h(t-tau) using Gauss-Chebyshev quadrature

  input:
    t:      time
    r:      parameter used in Bessel function argument
    fmax:   maximum frequency of Ricker wavelet
    mm:     size of array of Bessel function values
    func:   function used in numerical quadrature

  output:
    sum:    quadrature result
*/
{
  static int icall=0;
  int i, j, nn;
  float u, x, f;
  static double pi2;
  double pi, arg;
  static double *bes;

  void bessel_jn(double arg, int mm, double *bes);

  if (!icall) {
    icall = 1;
    bes = (double*)malloc((mm+1)*sizeof(double));
    pi  = 4.*atan(1.0);
    pi2 = pi/2.;
  }

  nn = t*r; /* number of quadrature points */
  nn /= 2;
  if (nn < 50) nn = 50;

  for (i=0; i<=mm; i++) sum[i] = 0.;

  for (i=0; i<nn; i++) {
    x = cos(pi2/nn*(2*i+1));

/* map interval [-1,1] to [0,t] */
    u = 0.5*t*(x+1.); 
    arg = u*r;

    f = (*func)(t-u,fmax) * pi/nn*sqrt(1.-x*x);

    bessel_jn(arg,mm,bes);

    for (j=1; j<=mm; j+=2) sum[j] += bes[j] * f;
  }

  /* normalize by interval length */
  for (j=1; j<=mm; j+=2) sum[j] *= 0.5 * t;

} /* end function intgr */
/**********************************************************************/
void intgru(float t, float r, float dt, int nt, int mm, float *wave, double* sum)
/*
  integrate J_n(tR) h(t-tau) using Gauss-Chebyshev quadrature

  input:
    t:      time
    r:      parameter used in Bessel function argument
    dt:     sample increment of user wavelet
    nt:     number of user wavelet samples 
    mm:     size of array of Bessel function values
    wave:   array of user wavelet values

  output:
    sum:    quadrature result
*/
{
  static int icall=0, nnmax=10000;
  int i, j, nn;
  float f;
  static float *x, *u, *xout, *yout;
  static double pi2;
  double pi, arg;
  static double *bes;

  void bessel_jn(double arg, int mm, double *bes);

  if (!icall) {
    icall = 1;
    bes = (double*)malloc((mm+1)*sizeof(double));
    x    = ealloc1float(nnmax);
    u    = ealloc1float(nnmax);
    xout = ealloc1float(nnmax);
    yout = ealloc1float(nnmax);
    pi  = 4.*atan(1.0);
    pi2 = pi/2.;
  }

  nn = t*r; /* number of quadrature points */
  nn /= 2;
  if (nn < 50) nn = 50;

  while (nn > nnmax) {
    fprintf(stderr,"nnmax=%d\n",nnmax);
    nnmax += 10000;
    if (nnmax > 100000) 
      throw cseis_geolib::csException("too many ( >100000 ) quadrature points in intgru!\n");
    x    = erealloc1float(x,nnmax);
    u    = erealloc1float(u,nnmax);
    xout = erealloc1float(xout,nnmax);
    yout = erealloc1float(yout,nnmax);
  }

  for (i=0; i<=mm; i++) sum[i] = 0.;
  for (i=0; i<nn; i++) {
    x[i] = cos(pi2/nn*(2*i+1));
    /* map interval [-1,1] to [0,t] */
    u[i] = 0.5*t*(x[i]+1.); 
    xout[i] = t-u[i];
  }

  /* sinc interpolation of wavelet to quadrature points */
  ints8r(nt,dt,0.,wave,wave[0],wave[nt-1],nn,xout,yout);
  for (i=0; i<nn; i++) {
    arg = u[i]*r;

    f = yout[i] * pi/nn*sqrt(1.-x[i]*x[i]);

    bessel_jn(arg,mm,bes);

    for (j=1; j<=mm; j+=2) sum[j] += bes[j] * f;
  }

  /* normalize by interval length */
  for (j=1; j<=mm; j+=2) sum[j] *= 0.5 * t;

} /* end function intgr */
/**********************************************************************/
float fwave(float t, float fmax)
/*
!
! --- Ricker wavelet ---
!
! t    : time
! fmax : maximum frequency of wavelet spectrum
!
! returns time sample at time t of a Ricker-wavelet with 
! maximum frequency fmax
!
*/
{
  static int ienter=0;
  static float pi, pi2, agauss, tcut, s, res;

  if (ienter == 0) {
    ienter = 1;
    pi = 4.*atan(1.);
    pi2 = 2.*pi;
    agauss = 0.5*fmax;
    tcut = 1.5/agauss;
  }

  s = (t-tcut) * agauss;

  if (fabs(s) < 4.) {
    res = exp(-2.*s*s) * cos(pi2*s);
  } else {
    res = 0.;
  }

  return res;
} /* end function fwave */
/**********************************************************************/
void damp(int nbw, float abso, float *bw)
{
/* 
    compute weights for aborbing boundary region

    input:
      nbw:    number of gridpoints of boundary width
      abso:   strength of absorption

    output:
      bw:     array of absorption weights

*/
  int i;
  float pi, delta;

  pi = 4. * atan(1.);
  delta = pi / nbw;
  
  for (i=0; i<nbw; i++) {
    bw[i] = 1.0 - abso * (1.0 + cos(i*delta)) * 0.5;
  }
 
  return;
}
/**********************************************************************/
void tbc2d(float *a, float *wbx, float *wbz, int nwbx, int nwbz,
           int nx, int nz)
/*
   2D absorbing boundary condition

   input:
   a:     array to be treated at boundaries
   wbx:   array of weights for absorbing boundaries in horizontal direction
   wbz:   array of weights for absorbing boundaries in vertical direction
   nbwx:  number of grid points of absorbing zone in x-direction
   nbwz:  number of grid points of absorbing zone in z-direction
   nx:    number of grid points of computational area in x-direction
   nz:    number of grid points of computational area in z-direction

   output:
   a:    array after absorbing boundaries treatment

   Note: The top absorbing region is placed below the bottom absorbing
         region. This is possible due to the cyclic nature of the
         discrete Fourier transform.

     (free) surface
   +----------------+
   |*              *|     * - absorbing region
   |*              *|
   |*              *|
   |*              *|
   |****************|
   |****************|
   +----------------+
         bottom
*/
{
  static int *ist;
  static int ienter = 0;
  int i, k;
  int ifl, ihv, istart, is, iz;
  int ind1, ind2, nwb, ng1, ng2, nz0;   
  float wb;
  
  if (ienter == 0) {
    ienter = 1;
    ist = ealloc1int(nz);
    ng1 = nwbx;
    ng2 = nwbz;
    ihv = 0;
    if (ng1 > ng2) {
      ng1 = nwbz;
      ng2 = nwbx;
      ihv = 1;
    }
    ifl = -ng2;
    istart = nwbx;
    iz = 0;
    if (ihv == 1) ist[0] = nwbx;
    for (i=0; i<ng2; i++) {
      if (ifl >= 0) {
	ifl = ifl-ng2;
	iz++;
	if (ihv == 1) ist[iz] = nwbx-i;
	istart--;
      }
      ifl = ifl + ng1;
      if (ihv == 0) ist[i] = istart;
      /* printf("i= %d ist= %d\n",i,ist[i]); */
    }
  }
/*
   +------------------------+
   |   Treatment of sides   |
   +------------------------+
*/
  nz0 = nz - 2*nwbz;
  for (k=0; k<nz; k++) {
    nwb = nwbx;
    if (k >= nz0+nwbz) nwb = ist[nz-k-1];
    if (k >= nz0 && k < nz0+nwbz) nwb = ist[k-nz0];

    for (i=0; i<nwb; i++) {
      wb = wbx[i];
      ind1 = i + k*nx;
      a[ind1] *= wb;
      ind2 = (nx-i-1) + k*nx;
      a[ind2] *= wb;
    }
  }

/*
   +-------------------------+
   |   Treatment of bottom   |
   +-------------------------+
*/
  for (k=0; k<nwbz; k++) {
    is = ist[k];
    wb = wbz[nwbz-k-1];

    for (i=is; i<nx-is; i++) {
      ind1 = i + (k+nz0)*nx;
      ind2 = i + (nz-k-1)*nx;
      a[ind1] *= wb;
      a[ind2] *= wb;
    }
  }
}

} // END namespace
