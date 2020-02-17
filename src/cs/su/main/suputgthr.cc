/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUPUTGTHR: $Revision: 1.5 $ ; $Date: 2011/11/17 00:03:38 $	*/

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
extern "C" {
  #include <unistd.h>

}

/*********************** self documentation **********************/
std::string sdoc_suputgthr =
" 									"
" SUPUTGTHR - split the stdout flow to gathers on the bases of given	"
" 		key parameter. 						"
" 									"
"	suputgthr <stdin   dir= [Optional parameters]			"
" 									"
" Required parameters:							"
" dir=		Name of directory where to put the gathers		"
" Optional parameters: 							"
" key=ep		header key word to watch   			"
" suffix=\".hsu\"	extension of the output files			"
" verbose=0		verbose = 1 echos information			"
" numlength=7		Length of numeric part of filename		"
" 									"
" Notes: 			    					"
" The name of the file is constructed from the key parameter. Traces	"
" are put into a temporary disk file, and renamed when key parameter	"
" changes in the input flow to \"key.suffix\". The result is that the	"
" directory \"dir\" contains separate files by \"key\" ensemble. 	"	
" 									"
" Header field modified:  ntr  to be the number of traces in a given 	"
" ensemble.								"
" Related programs: sugetgthr, susplit 					"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suputgthr {


/* 
 * Credits: Balazs Nemeth, Potash Corporation, Saskatoon Saskatchewan
 * given to CWP in 2008
 * Note:
 *	The "valxxx" subroutines are in su/lib/valpkge.c.  In particular,
 *	"valcmp" shares the annoying attribute of "strcmp" that
 *		if (valcmp(type, val, valnew) {
 *			...
 *		}
 *	will be performed when val and valnew are different.
 *
 */

/**************** end self doc ***********************************/


segy intrace;
segy tmptr;

void* main_suputgthr( void* args )
{
	cwp_String key;	/* header key word from segy.h		*/
	cwp_String type;/* ... its type				*/
	int index;	/* ... its index			*/
	int nsegy;	/* number of bytes in the segy		*/
	Value val;	/* value of key in current gather	*/
	Value valnew;	/* value of key in trace being treated	*/
	int verbose;	/* verbose flag				*/
	int val_i;	/* key value as an integer		*/

	int ntr=0;	/* count of number of traces in an ensemble */
	int numlength;	/* length of numerical part of filenames */
	
	FILE *tmpfp=NULL;		/* file pointer			*/
	FILE *outfp=NULL;	/* file pointer			*/
	cwp_String dir;		/* directory name		*/
	cwp_String suffix;	/* suffix of output files	*/

	char directory[BUFSIZ];		/* storage for directory name	*/
 	char tempfilename[BUFSIZ];	/* 	...temporary filename	*/
	char middle[BUFSIZ];		/*      ...middle of final filename */
	char stem[BUFSIZ];		/*      ...stem of final filename */
 	char format[BUFSIZ];		/* 	...format of numeric part */
	char outfile[BUFSIZ];		/*      ...final filename	*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suputgthr );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */



	/* Get parameters */
	CSMUSTGETPARSTRING("dir", &dir);
	if (!parObj.getparint   ("verbose", &verbose))	 	verbose = 0;
	if (!parObj.getparstring("key", &key))		 	key = "ep";
	if (!parObj.getparstring("suffix", &suffix))		suffix = ".hsu";
	if (!parObj.getparint("numlength", &numlength)) 	numlength=7;
        parObj.checkpars();

	/* Initialize */
	/* tempfilename = ealloc1(strlen(tmplt)+3,sizeof(char)); */
	type = hdtype(key);
	index = getindex(key);

	/* Set up for first trace (must compare new key field each time) */
	nsegy = cs2su->getTrace(&intrace);
	
	/* Create temporary filename for the first gather */
	strcpy(directory, dir);
	strcpy(tempfilename, temporary_filename(directory));
	if((tmpfp = fopen(tempfilename,"w+")) == NULL)
		throw cseis_geolib::csException("Cannot open file %s\n",tempfilename);

	if (verbose) warn(" temporary filename = %s", tempfilename);

	/* get key header value */
	gethval(&intrace, index, &val);
	
	ntr=0;
	do {

		++ntr;
		/* Get header value */				 
		gethval(&intrace, index, &valnew);

		/* compare past header value to current value */
		if (valcmp(type, val, valnew) || !nsegy) {
		 /* Either val and valnew differ, indicating a  */
		 /* new gather or nsegy is zero, indicating the */
		 /* end of the traces.			  */
			
			/* capture key field value */
			/* and build output filename */
			val_i=vtoi(type,val);

			/* zero out name parts */
			strcpy(outfile,"");
			strcpy(middle,"");
			strcpy(stem,"");

			/* build output name parts */
			strcat(middle,directory);
			strcat(middle,"/");
			strcat(stem,key);
			strcat(stem,"=");

			/* format for numeric part */
			(void)sprintf(format, "%%s%%s%%0%dd%%s",numlength);
		
			/* build name of output file */
			(void)sprintf(outfile, format,middle,stem,
					val_i, suffix);

			if (verbose) warn(" outfile = %s", outfile);

			/* rewind file */
			rewind(tmpfp);

			/* open the new file */
			if((outfp = fopen(outfile,"w+")) == NULL)
				throw cseis_geolib::csException("Cannot open file %s\n",outfile);
			/* loop over traces setting ntr field */
			/*  and write to finalfile */	
			while(fgettr(tmpfp,&tmptr))  {
				tmptr.ntr = ntr;
				fputtr(outfp,&tmptr);
			}
			/* Close files */
			fclose(tmpfp);
			fclose(outfp);
			remove(tempfilename);

			if (verbose) warn("val= %i", val_i);

			/* Set up for next gather */
			/* create new tempfname first */
			strcpy(tempfilename, temporary_filename(directory));
			
			/* open filename */
			if((tmpfp = fopen(tempfilename,"w+")) == NULL)
				throw cseis_geolib::csException("Cannot open file %s\n",tempfilename);
			val = valnew;
			ntr=0;
		}
		fputtr(tmpfp,&intrace);
	} while(cs2su->getTrace(&intrace));
	
	/* Close file */
	rewind(tmpfp);
	val_i=vtoi(type,val);
	
	/* Move last gather into new location */
	/* capture key field value */
	/* and build output filename */

	/* null out name parts */
	strcpy(outfile,"");
	strcpy(middle,"");
	strcpy(stem,"");

	/* build name parts */
	strcat(middle,directory);
	strcat(middle,"/");
	strcat(stem,key);
	strcat(stem,"=");
	
	/* write format of numeric part of output filename */
	(void)sprintf(format, "%%s%%s%%0%dd%%s",numlength);

	/* write output filename */
	(void)sprintf(outfile, format,middle,stem,
				val_i, suffix);

	/* open the output file */
	if((outfp = fopen(outfile,"w+")) == NULL)
			throw cseis_geolib::csException("Cannot open file %s\n",outfile);
	/* loop over traces setting ntr field */
	while(fgettr(tmpfp,&tmptr))  {
		tmptr.ntr = ntr;
		fputtr(outfp,&tmptr);
	}
	/* Close file */
	fclose(tmpfp);
	fclose(outfp);
	remove(tempfilename);

	if (verbose) warn(" outfile = %s", outfile);
	if (verbose) warn("val= %i",val_i);

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
