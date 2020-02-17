/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSUGetPars.h"
#include <cstdio>

using namespace cseis_su;
using namespace std;

csSUGetPars::csSUGetPars() {
  myXArgc = 0;
  myXArgv = NULL;
  myPar_count = 0;
  myParcheck = 0;
  myPar_names = new char*[cseis_su::PAR_NAMES_MAX];
  for( int i = 0; i < cseis_su::PAR_NAMES_MAX; i++ ) {
    myPar_names[i] = NULL;
  }

  myArgtbl = NULL;	/* parameter table		*/
  myNargs = 0;		/* number of args that parse	*/
  myTabled = FALSE;	        /* true when parameters tabled 	*/
  myTargc = 0;		/* total number of args		*/
  myTargv = NULL;		/* pointer to arg strings	*/
  myArgstr = NULL;		/* storage for command line	*/

}
csSUGetPars::~csSUGetPars() {
  if( myPar_names != NULL ) {
    for( int i = 0; i < myPar_count; i++ ) {
      if( myPar_names[i] != NULL ) {
        delete [] myPar_names[i];
      }
    }
    delete [] myPar_names;
    myPar_names = NULL;
  }
  if( myArgtbl != NULL ) {
    free(myArgtbl);
    myArgtbl = NULL;
  }
  if( myArgstr != NULL ) {
    free(myArgstr);
    myArgstr = NULL;
  }
}

void csSUGetPars::addPar( char const* text ) {
  int length = strlen(text);
  myPar_names[myPar_count] = new char[length+1];
  memcpy( myPar_names[myPar_count], text, length );
  myPar_names[myPar_count][length] = '\0';
  myPar_count += 1;
}

void csSUGetPars::initargs( int argc, char **argv )
{
  addPar( "par" );
  addPar( "lheader" );
  
  myXArgc = argc;
  myXArgv = argv;
  if(myTabled==TRUE){
    free(myArgstr);
    free(myTargv);
    free(myArgtbl);
  }
  myTabled =  FALSE;
}

int csSUGetPars::getparint (char *name, unsigned int *ptr)
{
  int valueInt;
  int result = getnpar(0,name,"i",&valueInt);
  *ptr = (unsigned int)valueInt;
  return result;
}
int csSUGetPars::getparint (char *name, int *ptr)
{
	return getnpar(0,name,"i",ptr);
}
int csSUGetPars::getparuint (char *name, unsigned int *ptr)
{
	return getnpar(0,name,"p",ptr);
}
int csSUGetPars::getparshort (char *name, short *ptr)
{
	return getnpar(0,name,"h",ptr);
}
int csSUGetPars::getparushort (char *name, unsigned short *ptr)
{
	return getnpar(0,name,"u",ptr);
}
int csSUGetPars::getparlong (char *name, long *ptr)
{
	return getnpar(0,name,"l",ptr);
}
int csSUGetPars::getparulong (char *name, unsigned long *ptr)
{
	return getnpar(0,name,"v",ptr);
}
int csSUGetPars::getparfloat (char *name, float *ptr)
{
	return getnpar(0,name,"f",ptr);
}
int csSUGetPars::getpardouble (char *name, double *ptr)
{
	return getnpar(0,name,"d",ptr);
}
int csSUGetPars::getparstring (char *name, char **ptr)
{
	return getnpar(0,name,"s",ptr);
}
int csSUGetPars::getparstringarray (char *name, char **ptr)
{
	return getnpar(0,name,"a",ptr);
}
int csSUGetPars::getpar(char *name, char *type, void *ptr)
{
	return getnpar(0,name,type,ptr);
}

int csSUGetPars::getnparint (int n, char *name, int *ptr)
{
	return getnpar(n,name,"i",ptr);
}
int csSUGetPars::getnparuint (int n, char *name, unsigned int *ptr)
{
	return getnpar(n,name,"p",ptr);
}
int csSUGetPars::getnparshort (int n, char *name, short *ptr)
{
	return getnpar(n,name,"h",ptr);
}
int csSUGetPars::getnparushort (int n, char *name, unsigned short *ptr)
{
	return getnpar(n,name,"u",ptr);
}
int csSUGetPars::getnparlong (int n, char *name, long *ptr)
{
	return getnpar(n,name,"l",ptr);
}
int csSUGetPars::getnparulong (int n, char *name, unsigned long *ptr)
{
	return getnpar(n,name,"v",ptr);
}
int csSUGetPars::getnparfloat (int n, char *name, float *ptr)
{
	return getnpar(n,name,"f",ptr);
}
int csSUGetPars::getnpardouble (int n, char *name, double *ptr)
{
	return getnpar(n,name,"d",ptr);
}
int csSUGetPars::getnparstring (int n, char *name, char **ptr)
{
	return getnpar(n,name,"s",ptr);
}
int csSUGetPars::getnparstringarray (int n, char *name, char **ptr)
{
	return getnpar(n,name,"a",ptr);
}
int csSUGetPars::getnpar (int n, char *name, char *type, void *ptr)
{
	int i;			/* index of name in symbol table	*/
	int j;		  /* index for myPar_names[]		*/
	int nval;		/* number of parameter values found	*/
	char *aval;		/* ascii field of symbol		*/

/*--------------------------------------------------------------------*\
   getpar gets called in loops reading traces in some programs.  So
   check for having seen this name before. Also make sure we don't
   walk off the end of the table.
\*--------------------------------------------------------------------*/

	if( myParcheck && strcmp( "lheader" ,name ) ){
	   fprintf( stderr ,"csSUGetPars::getpar() call after checkpars(): %s\n" ,name );
	}

	for( j=0; j<myPar_count; j++ ){
	   if( !strcmp( myPar_names[j] ,name ) ){
		break;
	   }
	}

	if( j >= myPar_count && myPar_count < cseis_su::PAR_NAMES_MAX ){
          addPar( name );
	}

	if(  myPar_count == cseis_su::PAR_NAMES_MAX ){
	   fprintf( stderr, " %s exceeded PAR_NAMES_MAX %d \n" ,myXArgv[0] ,cseis_su::PAR_NAMES_MAX );
	}

	if (myXArgc == 1) return 0;
	if (!myTabled) getparinit();/* Tabulate command line and parfile */
	i = getparindex(n,name);/* Get parameter index */
	if (i < 0) return 0;	/* Not there */

	if (0 == ptr) {
	   err("%s: getnpar called with 0 pointer, type = %s", __FILE__,type);
	}
	  

	/*
	 * handle string type as a special case, since a string
	 * may contain commas.
	 */
	if (type[0]=='s') {
		*((char**)ptr) = myArgtbl[i].asciival;
		return 1;
	}

	/* convert vector of ascii values to numeric values */
	for (nval=0,aval=myArgtbl[i].asciival; *aval; nval++) {
		switch (type[0]) {
			case 'i':
				*(int*)ptr = eatoi(aval);
				ptr = (int*)ptr+1;
				break;
			case 'p':
				*(unsigned int*)ptr = eatop(aval);
				ptr = (unsigned int*)ptr+1;
				break;
			case 'h':
				*(short*)ptr = eatoh(aval);
				ptr = (short*)ptr+1;
				break;
			case 'u':
				*(unsigned short*)ptr = eatou(aval);
				ptr = (unsigned short*)ptr+1;
				break;
			case 'l':
				*(long*)ptr = eatol(aval);
				ptr = (long*)ptr+1;
				break;
			case 'v':
				*(unsigned long*)ptr = eatov(aval);
				ptr = (unsigned long*)ptr+1;
				break;
			case 'f':
				*(float*)ptr = eatof(aval);
				ptr = (float*)ptr+1;
				break;
			case 'd':
				*(double*)ptr = eatod(aval);
				ptr = (double*)ptr+1;
				break;
			case 'a':
				{ char *tmpstr="";
                                  tmpstr =(char*) ealloc1(strlen(aval)+1,1);

				   strchop(aval,tmpstr);
				   *(char**)ptr = tmpstr;
				   ptr=(char **)ptr + 1;
				}
				   break;
			default:
				err("%s: invalid parameter type = %s",
					__FILE__,type);
		}
		while (*aval++ != ',') {
			if (!*aval) break;
		}
	}
	return nval;
}

void csSUGetPars::checkpars(){

   int i;
   int j;
   char buf[256];

   for( j=1; j<myXArgc; j++){

	for( i=0; i<myPar_count; i++ ){
	 sprintf( buf ,"%s=" ,myPar_names[i] );

	 if( !strncmp( buf ,myXArgv[j] ,strlen(buf) ) ){
	    break;
	 }
	}
	if( i == myPar_count && strchr( myXArgv[j] ,'=' ) ){
	 fprintf( stderr ,"Unknown %s argument %s\n" ,myXArgv[0] ,myXArgv[j] );
	}

   }

   myParcheck = 1;
}

/* Promax compatible version of getnpar */
void csSUGetPars::getparPromax(char *name, char *type, void *ptr)
{
	(void) getnpar(0,name,type,ptr);
	return;
}

/* return number of occurrences of parameter name */
int csSUGetPars::countparname (char *name)
{
	int i,nname;

	if (myXArgc == 1) return 0;
	if (!myTabled) getparinit();
	for (i=0,nname=0; i<myNargs; ++i)
		if (!strcmp(name,myArgtbl[i].name)) ++nname;
	return nname;
}

/* return number of values in n'th occurrence of parameter name */
int csSUGetPars::countnparval (int n, char *name)
{
	int i;

	if (myXArgc == 1) return 0;
	if (!myTabled) getparinit();
	i = getparindex(n,name);
	if (i>=0)
		return ccount(',',myArgtbl[i].asciival) + 1;
	else
		return 0;
}

/* return number of values in last occurrence of parameter name */
int csSUGetPars::countparval (char *name)
{
	return countnparval(0,name);
}



/*
 * Return the index of the n'th occurrence of a parameter name,
 * except if n==0, return the index of the last occurrence.
 * Return -1 if the specified occurrence does not exist.
 */
int csSUGetPars::getparindex (int n, char *name)
{
	int i;
	if (n==0) {
		for (i=myNargs-1; i>=0; --i)
			if (!strcmp(name,myArgtbl[i].name)) break;
		return i;
	} else {
		for (i=0; i<myNargs; ++i)
			if (!strcmp(name,myArgtbl[i].name))
				if (--n==0) break;
		if (i<myNargs)
			return i;
		else
			return -1;
	}
}

/* Initialize getpar */

void csSUGetPars::getparinit (void)
{
	static char *pfname;	/* name of parameter file		*/
	FILE *pffd=NULL;	/* file id of parameter file		*/
	size_t pflen;		/* length of parameter file in bytes	*/
	int parfile;		/* parfile existence flag		*/
	int myArgstrlen=0;
	char *pmyArgstr;		/* storage for parameter file args	*/
	size_t nread=0;		/* bytes fread				*/
	int i, j;		/* counters				*/
	int start = TRUE;
	int debug = FALSE;
	int quote = FALSE;


	myTabled = TRUE;		/* remember table is built		*/


	/* Check if myXArgc was initiated */

	if(!myXArgc)
		err("%s: myXArgc=%d -- not initiated in main", __FILE__, myXArgc);

	/* Space needed for command lines */

	for (i = 1, myArgstrlen = 0; i < myXArgc; i++) {
		myArgstrlen += strlen(myXArgv[i]) + 1;
	}

	/* Get parfile name if there is one */

	if ((pfname = getpfname())) {
		parfile = TRUE;
	} else {
		parfile = FALSE;
	}

	if (parfile) {
	 	pffd = efopen(pfname, "r");

		/* Get the length */
		efseek(pffd, 0, SEEK_END);

		pflen = (off_t) eftello(pffd);

		rewind(pffd);
		myArgstrlen += pflen;
	} else {
		pflen = 0;
	}

/*--------------------------------------------------------------------*\
   Allocate space for command line and parameter file. The pointer
   table could be as large as the string buffer, but no larger.

   The parser logic has been completely rewritten to prevent bad
   input from crashing the program.

   Reginald H. Beardsley			    rhb@acm.org
\*--------------------------------------------------------------------*/

	myArgstr = (char *) ealloc1(myArgstrlen+1, 1);
	myTargv = (char **) ealloc1((myArgstrlen+1)/4,sizeof(char*));

	if (parfile) {
		/* Read the parfile */

		nread = efread(myArgstr, 1, pflen, pffd);
  		if (nread != pflen) {
  	 	    err("%s: fread only %d bytes out of %d from %s",
  					__FILE__,  nread, pflen, pfname);
		}
		efclose(pffd);


	}


	/* force input to valid 7 bit ASCII */

	for( i=0; i<(int)nread; i++ ){
	    myArgstr[i] &= 0x7F;
	}

	/* tokenize the input */

	j = 0;

	for( i=0; i<(int)nread; i++ ){

	    /* look for start of token */

	    if( start ){

 /* getpars.c:475: warning: subscript has type `char' */
		if( isgraph( (int)myArgstr[i] ) ){
		    myTargv[j] = &(myArgstr[i]);
		    start = !start;
		    j++;

		}else{
		    myArgstr[i] = 0;

		}

	    /* terminate token */

/* getpars.c:487: warning: subscript has type `char' */
	    }else if( !quote && isspace( (int)myArgstr[i] ) ){
		myArgstr[i] = 0;
		start = !start;

	    }

	    /* toggle quote semaphore */

	    if( myArgstr[i] == '\'' || myArgstr[i] == '\"' ){
		quote = !quote;

	    }

	}

	/* display all tokens */

	if( debug ){

	    i=0;
	    while( i < j && myTargv[i] != 0 ){
		if( strlen( myTargv[i] ) ){
		    fprintf( stderr ,"%d -> %s\n" ,i ,myTargv[i] );
		}
		i++;

	    }
	}

	/* discard non-parameter tokens */

	i=0;
	myTargc=0;
	while( i < j && myTargv[i] != 0 ){
	    if( strchr( myTargv[i] ,'=' ) ){
		myTargv[myTargc] = myTargv[i];
		myTargc++;
	    }
	    i++;
	}

	/* Copy command line arguments */

	for (j = 1, pmyArgstr = myArgstr + pflen + 1; j < myXArgc; j++) {
		strcpy(pmyArgstr,myXArgv[j]);
		myTargv[myTargc++] = pmyArgstr;
		pmyArgstr += strlen(myXArgv[j]) + 1;
	}

	/* Allocate space for the pointer table */

	myArgtbl = (cseis_su::pointer_table*) ealloc1(myTargc, sizeof(cseis_su::pointer_table));

	/* Tabulate myTargv */

	tabulate(myTargc, myTargv);

	return;
}

char* csSUGetPars::getpfname (void)
{
	int i;
	size_t pfnamelen;

	pfnamelen = strlen("par=");
	for (i = myXArgc-1 ; i > 0 ; i--) {
		if(!strncmp("par=", myXArgv[i], pfnamelen)
		    && strlen(myXArgv[i]) != pfnamelen) {
			return myXArgv[i] + pfnamelen;
		}
	}
	return NULL;
}

int iswhite(char letter) {
  return( ((letter) == ' ' || (letter) == '\t' || (letter) == '\n') );
}

size_t csSUGetPars::white2null (char *str, size_t len)
{
	int i;
	size_t count = 0;
	int inquote = FALSE;

	str[0] = '\0'; /* This line added by Dave Hale, 1/30/96. */
	for (i = 1; i < (int)len; i++) {
		if (str[i]=='"') inquote=(inquote==TRUE)?FALSE:TRUE;
		if (!inquote) {
			if (iswhite(str[i])) { /* Is this a new word ? */
				str[i] = '\0';
			} else if (!str[i-1]) { /* multiple whites */
				count++;
			}
		}
	}
	for (i = 1, inquote=FALSE; i < (int)len; i++) {
		if (str[i]=='"') inquote=(inquote==TRUE)?FALSE:TRUE;
		if (inquote) {
			if (str[i+1]!='"') {
				str[i] = str[i+1];
			} else {
				str[i] = '\0';
				str[i+1] = '\0';
				inquote = FALSE;
			}
		}
	}
	str[len] = '\0';
	return count;
}

/* Install symbol table */
void csSUGetPars::tabulate (size_t argc, char **argv)
{
	int i;
	char *eqptr;
	int debug=FALSE;

	for (i = 0, myNargs = 0 ; i < (int)argc; i++) {
		eqptr = strchr(argv[i], '=');
		if (eqptr) {
			myArgtbl[myNargs].name = argv[i];
			myArgtbl[myNargs].asciival = eqptr + 1;
			*eqptr = (char)0;

			/* Debugging dump */
			if( debug ){
				fprintf(stderr,
				"myArgtbl[%d]: name=%s asciival=%s\n",
				myNargs,myArgtbl[myNargs].name,myArgtbl[myNargs].asciival);

			}
			myNargs++;
		}
	}
	return;
}

/* Count characters in a string */
int csSUGetPars::ccount (char c, char *s)
{
	int i, count;
	for (i = 0, count = 0; s[i] != 0; i++)
		if(s[i] == c) count++;
	return count;
}
