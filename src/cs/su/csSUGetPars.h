/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/*****************************************************************************
Authors:
Rob Clayton & Jon Claerbout, Stanford University, 1979-1985
Shuki Ronen & Jack Cohen, Colorado School of Mines, 1985-1990
Dave Hale, Colorado School of Mines, 05/29/90
Credit to John E. Anderson for re-entrant initargs 03/03/94
*****************************************************************************/
/**************** end self doc ********************************/
#ifndef CS_SU_GET_PARS_H
#define CS_SU_GET_PARS_H

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "par.h"

/*--------------------------------------------------------------------*\
  These variables are used by checkpars() to warn of parameter typos.
  par= does not use getpar() so we need to store that as the first
  parameter name.  lheaders= is buried in fgettr() so we intialize
  that also
\*--------------------------------------------------------------------*/

#define CSMUSTGETPARINT(x,y)    if(!parObj.getparint(x,y))    { throw(cseis_geolib::csException("SU module %s: Must specify %s=\n",argv[0],x)); }
#define CSMUSTGETPARFLOAT(x,y)  if(!parObj.getparfloat(x,y))  { throw(cseis_geolib::csException("SU module %s: Must specify %s=\n",argv[0],x)); }
#define CSMUSTGETPARSTRING(x,y) if(!parObj.getparstring(x,y)) { throw(cseis_geolib::csException("SU module %s: Must specify %s=\n",argv[0],x)); }
#define CSMUSTGETPARDOUBLE(x,y) if(!parObj.getpardouble(x,y)) { throw(cseis_geolib::csException("SU module %s: Must specify %s=\n",argv[0],x)); }

namespace cseis_su {

  static int const PAR_NAMES_MAX = 512;

  /* par fprintf(stderr,"\n");su2cs->setEOF();ameter table */
  struct pointer_table {
    char *name;      /* external name of parameter */
    char *asciival;  /* ascii value of parameter   */
  };

class csSUGetPars {
 public:
  csSUGetPars();
  ~csSUGetPars();
  void initargs( int argc, char **argv );
  int getparint (char *name, unsigned int *ptr);
  int getparint (char *name, int *ptr);
  int getparuint (char *name, unsigned int *ptr);
  int getparshort (char *name, short *ptr);
  int getparushort (char *name, unsigned short *ptr);
  int getparlong (char *name, long *ptr);
  int getparulong (char *name, unsigned long *ptr);
  int getparfloat (char *name, float *ptr);
  int getpardouble (char *name, double *ptr);
  int getparstring (char *name, char **ptr);
  int getparstringarray (char *name, char **ptr);
  int getpar (char *name, char *type, void *ptr);
  int getnparint (int n, char *name, int *ptr);
  int getnparuint (int n, char *name, unsigned int *ptr);
  int getnparshort (int n, char *name, short *ptr);
  int getnparushort (int n, char *name, unsigned short *ptr);
  int getnparlong (int n, char *name, long *ptr);
  int getnparulong (int n, char *name, unsigned long *ptr);
  int getnparfloat (int n, char *name, float *ptr);
  int getnpardouble (int n, char *name, double *ptr);
  int getnparstring (int n, char *name, char **ptr);
  int getnparstringarray (int n, char *name, char **ptr);
  int getnpar (int n, char *name, char *type, void *ptr);
  void checkpars();
  void getparPromax(char *name, char *type, void *ptr);
  int countparname (char *name);
  int countnparval (int n, char *name);
  int countparval (char *name);
  int getparindex (int n, char *name);
  void getparinit (void);
  
  char *getpfname (void);
  void tabulate (size_t argc, char **argv);

  static size_t white2null (char *str, size_t len);
  static int ccount (char c, char *s);
  static int iswhite(char letter);

 private:
  void addPar( char const* text );

  char** myPar_names;
  int   myPar_count;
  int   myParcheck;
  int   myXArgc;
  char** myXArgv;


  cseis_su::pointer_table *myArgtbl; /* parameter table */
  int myNargs;          /* number of args that parse    */
  int myTabled;         /* true when parameters tabled  */
  size_t myTargc;       /* total number of args         */
  char **myTargv;       /* pointer to arg strings       */
  char *myArgstr;       /* storage for command line     */
};

} // end namespace

#endif
