/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUGETGTHR: $Revision: 1.5 $ ; $Date: 2011/11/17 00:03:38 $		*/

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
  #include <sys/stat.h> 

}
extern "C" {
  #include <sys/types.h>

}
extern "C" {
  #include <fcntl.h>

}
extern "C" {
  #include <dirent.h>

}

#include "su.h"
#include "segy.h"
#include "header.h"


#define CWP_O_LARGEFILE	0100000

/*********************** self documentation **********************/
std::string sdoc_sugetgthr =
" 									"
" SUGETGTHR - Gets su files from a directory and put them               "
"             throught the unix pipe. This creates continous data flow.	"
" 									"
"  sugetgthr  <stdin >sdout   						"
" 									"
" Required parameters:							"
" 									"
" dir=            Name of directory to fetch data from 			"
" 	          Every file in the directory is treated as an su file	"
" Optional parameters:							"
" verbose=0		=1 more chatty					"
" vt=0			=1 allows gathers with variable length traces	"
" 			no header checking is done!			"
" ns=			must be specified if vt=1; number of samples to read"
" 									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace sugetgthr {




segy tr;

void* main_sugetgthr( void* args )
{
	
	cwp_String dir="";	/* input directory containng the gathers */
	char *fname=NULL;	
	char *ffname=NULL;
	
	DIR *dp=NULL;
	struct dirent *d=NULL;
	struct stat __st;
	FILE *fp=NULL;
	int fd=0;
	
	int verbose;
	int vt=0;
	ssize_t nread;
	int ns=0;

	
	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_sugetgthr );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


        /* Get parameters */
        CSMUSTGETPARSTRING("dir", &dir);
       	if (!parObj.getparint   ("verbose", &verbose)) verbose = 0;
       	if (!parObj.getparint   ("vt", &vt)) vt = 0;
	if(vt)CSMUSTGETPARINT("ns",&ns); 
        parObj.checkpars();
	
	/* Open the directory */
	if ((dp = opendir(dir)) == NULL)
		throw cseis_geolib::csException(" %s directory not found\n",dir);
	
	/* For each file in directory */
	while (( d = readdir(dp)) !=NULL) {
		
		fname = ealloc1(strlen(d->d_name)+1,sizeof(char));
		strcpy(fname,d->d_name);
		
		/* Skip . and .. directory entries */
		if(strcmp(fname,".") && strcmp(fname,"..")) {		
			ffname = ealloc1(strlen(d->d_name)+strlen(dir)+2,sizeof(char));
			
			/* Create full filename */
			sprintf(ffname, "%s/%s",dir,fname);
			if(verbose==1) warn("%s",ffname);
			
			/* get some info from the file */
			stat(ffname,&__st);
			if(__st.st_size > 0) {
			
				/* Open the file and read traces into stdout*/
 				if(vt) {
					fd = open(ffname,O_RDONLY|CWP_O_LARGEFILE);
				/*	nread=fread(&tr,(size_t) HDRBYTES,1,fp);  */
					nread=read(fd,&tr,(size_t) HDRBYTES); 
					memset((void *) &tr.data[tr.ns], (int) '\0' ,MAX(ns-tr.ns,0)*FSIZE);
				/*	nread+=fread(&tr.data[0],(size_t) tr.ns*FSIZE,1,fp); */
					nread+=read(fd,&tr.data[0],(size_t) tr.ns*FSIZE);
				} else {
					fp = efopen(ffname, "r");
					nread=fgettr(fp, &tr);
				}
				do {
					if(vt) { 
						tr.ns=ns;
						fwrite(&tr,ns*FSIZE+HDRBYTES,1,stdout);
					} else {
						su2cs->putTrace(&tr);
					}
					if(vt) {
						/* nread=fread(&tr,(size_t) HDRBYTES,1,fp); */
						nread=read(fd,&tr,(size_t) HDRBYTES);
						memset((void *) &tr.data[tr.ns], (int) '\0' ,MAX(ns-tr.ns,0)*FSIZE);
						/* nread+=fread(&tr.data[0],(size_t) tr.ns*FSIZE,1,fp); */
						nread+=read(fd,&tr.data[0],(size_t) tr.ns*FSIZE);
					} else {
						nread=fgettr(fp, &tr);
					}
				} while(nread);
				if(vt) close(fd);
				else efclose(fp);
			} else {
				warn(" File %s has zero size, skipped.\n",ffname);
			}
			free1(ffname);
		}
		free1(fname);
		
	}
	closedir(dp);
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
