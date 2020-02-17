/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUMIGSPLIT: $Revision: 1.9 $ ; $Date: 2011/11/16 22:14:43 $       */

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
extern "C" {
  #include <signal.h>

}

/*********************** self documentation ******************************/
std::string sdoc_sumigsplit =
"                                                                       "
" SUMIGSPLIT - Split-step depth migration for zero-offset data.         "
"                                                                       "
" sumigsplit <infile >outfile vfile= [optional parameters]              "
"                                                                       "
" Required Parameters:                                                  "
" nz=                   number of depth sapmles                         "
" dz=                   depth sampling interval                         "
" vfile=                name of file containing velocities              "
"                                                                       "
" Optional Parameters:                                                  "
" dt=from header(dt) or .004    time sampling interval                  "
" dx=from header(d2) or 1.0     midpoint sampling interval              "
" ft=0.0                        first time sample                       "
" fz=                           first depth sample                      "
"                                                                       "
" tmpdir=        if non-empty, use the value as a directory path        "
"                prefix for storing temporary files; else if the        "
"                the CWP_TMPDIR environment variable is set use         "
"                its value for the path; else use tmpfile()             "
"                                                                       "
"                                                                       "
" Notes:                                                                "
" The input velocity file \'vfile\' consists of C-style binary floats.  "
" The structure of this file is vfile[iz][ix]. Note that this means that"
" the x-direction is the fastest direction instead of z-direction! Such a"
" structure is more convenient for the downward continuation type       "
" migration algorithm than using z as fastest dimension as in other SU  "
" programs.                                                     	"
"                                                                       "
" Because most of the tools in the SU package (such as  unif2, unisam2, "
" and makevel) produce output with the structure vfile[ix][iz], you will"
" need to transpose the velocity files created by these programs. You may"
" use the SU program \'transp\' in SU to transpose such files into the  "
" required vfile[iz][ix] structure.                                     "
" (In C  v[iz][ix] denotes a v(x,z) array, whereas v[ix][iz]  		"
" denotes a v(z,x) array, the opposite of what Matlab and Fortran	"
" programmers may expect.)						" 
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sumigsplit {

  
/*
 * Credits: CWP Baoniu Han, July 21th, 1997
 *
 * Reference:
 *  Stoffa, P. L. and Fokkema, J. T. and Freire, R. M. and Kessinger, W. P.,
 *  1990, Split-step Fourier migration, Geophysics, 55, 410-421.
 *
 *
 * Trace header fields accessed: ns, dt, delrt, d2
 * Trace header fields modified: ns, dt, delrt
 */

/**************** end self doc *******************************************/

/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
static void closefiles(void);
segy tr;



void* main_sumigsplit( void* args )

{
	int nt;			/* number of time samples */
	int nz;			/* number of migrated depth samples */
	int nx;			/* number of midpoints 	*/
	int ik,iz,iw,ix,it;	/* loop counters 	*/
	int nxfft,ntfft;	/* fft size		*/
	int nk,nw;		/* number of wave numbers */	


	float dt,dz,tz;		/* time sampling interval 	*/
	float ft,fz;		/* first time sample		*/
	float dk,dw;		/* wave number sampling interval */
	float fk,fw;		/* first wave number 		*/
	float k,w;		/* wave number		*/
	float dx;		/* spatial sampling interval	*/
	float **p,**cresult;	/* input, output data		*/
	float vmin,v1;
	double kz1,kz2;
	float **v;
	double phase1;
	complex cshift1,cshift2;
	complex **cp,**cq,**cq1;	/* complex input,output		*/
	char *vfile="";		/* name of file containing velocities */
	FILE *vfp;
	int verbose=1;		/* flag for echoing info		*/
	char *tmpdir;		/* directory path for tmp files		*/
	cwp_Bool istmpdir=cwp_false;/* cwp_true for user-given path		*/

	/* hook up getpar to handle the parameters */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sumigsplit );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* get info from first trace */
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	nt = tr.ns;

	/* let user give dt and/or dx from command line */
	if (!parObj.getparfloat("dt", &dt)) {
		if (tr.dt) { /* is dt field set? */
			dt = ((double) tr.dt)/1000000.0;
		} else { /* dt not set, assume 4 ms */
			dt = 0.004;
			warn("tr.dt not set, assuming dt=0.004");
		}
	}
	if (!parObj.getparfloat("dx",&dx)) {
		if (tr.d2) { /* is d2 field set? */
			dx = tr.d2;
		} else {
			dx = 1.0;
			warn("tr.d2 not set, assuming d2=1.0");
		}
	}


	/* get optional parameters */
	if (!parObj.getparfloat("ft",&ft)) ft = 0.0;
        if (!parObj.getparint("nz",&nz)) throw cseis_geolib::csException("nz must be specified");
        if (!parObj.getparfloat("dz",&dz)) throw cseis_geolib::csException("dz must be specified");
        if (!parObj.getparfloat("fz",&fz)) fz = 0.0;
        if (!parObj.getparstring("vfile", &vfile)) throw cseis_geolib::csException("vfile must be specified");

	if (!parObj.getparint("verbose", &verbose)) verbose = 0;

	/* Look for user-supplied tmpdir */
	if (!parObj.getparstring("tmpdir",&tmpdir) &&
	    !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
	if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
		throw cseis_geolib::csException("you can't write in %s (or it doesn't exist)", tmpdir);
	/* store traces and headers in tempfiles while getting a count */
	if (STREQ(tmpdir,"")) {
		tracefp = etmpfile();
		headerfp = etmpfile();
		if (verbose) warn("using tmpfile() call");
	} else { /* user-supplied tmpdir */
		char directory[BUFSIZ];
		strcpy(directory, tmpdir);
		strcpy(tracefile, temporary_filename(directory));
		strcpy(headerfile, temporary_filename(directory));
		/* Trap signals so can remove temp files */
		signal(SIGINT,  (void (*) (int)) closefiles);
		signal(SIGQUIT, (void (*) (int)) closefiles);
		signal(SIGHUP,  (void (*) (int)) closefiles);
		signal(SIGTERM, (void (*) (int)) closefiles);
		tracefp = efopen(tracefile, "w+");
		headerfp = efopen(headerfile, "w+");
      		istmpdir=cwp_true;		
		if (verbose) warn("putting temporary files in %s", directory);
	}

        parObj.checkpars();

	nx = 0;
	do {
		 ++nx;
		efwrite(&tr,HDRBYTES,1,headerfp);
		efwrite(tr.data, FSIZE, nt, tracefp);
	} while (cs2su->getTrace(&tr));
	erewind(tracefp);
	erewind(headerfp);
	
	/* determine frequency sampling (for real to complex FFT) */
	ntfft = npfar(nt);
	nw = ntfft/2+1;
	dw = 2.0*PI/(ntfft*dt);
	fw = 0;
	
	/* determine wavenumber sampling (for complex to complex FFT) */
	nxfft = npfa(nx);
	nk = nxfft;
	dk = 2.0*PI/(nxfft*dx);
	fk = -PI/dx;


	/* allocate space */
	p = alloc2float(ntfft,nx);
	cp = alloc2complex(nw,nx);
	cq = alloc2complex(nw,nk);
	cq1 = alloc2complex(nw,nk);
	cresult = alloc2float(nz,nx);
	v=alloc2float(nx,nz);

	/* load traces into the zero-offset array and close tmpfile */
	for(ix=0;ix<nx;ix++)
	efread(p[ix], FSIZE, nt, tracefp);
	efclose(tracefp);

	/*load velicoty file*/
	vfp=efopen(vfile,"r");	
	efread(v[0],FSIZE,nz*nx,vfp);
	efclose(vfp);			
	
	/* pad with zeros and Fourier transform x to k */
	for (ix=0; ix<nx; ix++)
		for (it=nt; it<ntfft; it++)
			p[ix][it] = 0.0;

	pfa2rc(1,1,ntfft,nx,p[0],cp[0]);

	/* loops over depth */
	for(iz=0,tz=fz;iz<nz;++iz,tz+=dz){
		for(ix=0;ix<nx;ix++){
			cresult[ix][iz] =0.0;
			for(iw=0;iw<nw;iw++)
			cresult[ix][iz]+=cp[ix][iw].r/ntfft;

			for(iw=0;iw<nw;iw++)
			cp[ix][iw].r-=cresult[ix][iz]/ntfft;


		}
		
		vmin=0;

		for(ix=0;ix<nx;ix++){ 
		vmin+=1.0/v[iz][ix]/nx;
		}

		vmin=1.0/vmin*0.5;

		for (ik=0; ik<nx; ++ik)
		        for (iw=0; iw<nw; ++iw) 
			        cq[ik][iw] = ik%2 ? cneg(cp[ik][iw]) : cp[ik][iw];

		for (ik=nx; ik<nk; ++ik)
			for (iw=0; iw<nw; ++iw)
				cq[ik][iw] = cmplx(0.0,0.0);


		/* FFT to W-K domain */

		pfa2cc(-1,2,nw,nk,cq[0]);

		v1=vmin;

		for(ik=0,k=fk;ik<nk;++ik,k+=dk)
			for(iw=0,w=fw;iw<nw;++iw,w+=dw){
				if(w==0.0)w=1.0e-10/dt;
				kz1=1.0-pow(v1*k/w,2.0);
				if(kz1>0){
				phase1 = -w*sqrt(kz1)*dz/v1; 
				cshift1 = cmplx(cos(phase1), sin(phase1));
				cq1[ik][iw] = cmul(cq[ik][iw],cshift1);
				}
				else{
				cshift1=cmplx(0.0,0.0);
				cq1[ik][iw] = cmul(cq[ik][iw],cshift1);
				}

			}


		pfa2cc(1,2,nw,nk,cq1[0]);

		for(ix=0;ix<nx;++ix)
			for(iw=0,w=fw;iw<nw;w+=dw,++iw){
				cq1[ix][iw] = crmul( cq1[ix][iw], 1.0/nxfft);
				cq1[ix][iw] =ix%2 ? cneg(cq1[ix][iw]) : cq1[ix][iw];
				kz2=(1.0/v1-2.0/v[iz][ix])*w*dz;
				cshift2=cmplx(cos(kz2),sin(kz2));
				cp[ix][iw]=cmul(cq1[ix][iw],cshift2);
			}



		}
	/* restore header fields and write output */
	for(ix=0; ix<nx; ix++){
		efread(&tr,HDRBYTES,1,headerfp);
		tr.ns = nz ;
		tr.d1 = dz ;
		memcpy( (void *) tr.data, (const void *) cresult[ix],nz*FSIZE);
		su2cs->putTrace(&tr);
	}
	
	/* Clean up */
	efclose(headerfp);
	if (istmpdir) eremove(headerfile);
	if (istmpdir) eremove(tracefile);

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


static void closefiles(void)
{
        efclose(headerfp);
        efclose(tracefp);
        eremove(headerfile);
        eremove(tracefile);
        exit(EXIT_FAILURE);
}


} // END namespace
