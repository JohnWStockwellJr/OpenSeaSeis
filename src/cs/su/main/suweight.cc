/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUWEIGHT: $Revision: 1.8 $ ; $Date: 2011/11/16 17:23:05 $   */

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

/*************************** self documentation ******************************/
std::string sdoc_suweight =
"									"
" SUWEIGHT - weight traces by header parameter, such as offset		"
"									"
"   suweight < stdin > stdout [optional parameters]			"
"									"
" Required Parameters:					   		"
"   <none>								"
"									"
" Optional parameters:					   		"
" key=offset	keyword of header field to weight traces by 		"
" a=1.0		constant weighting parameter (see notes below)		"
" b=.0005	variable weighting parameter (see notes below)		"
"									"
"... or use values of a header field for the weighting ...		"
"									"
" key2=		keyword of header field to draw weights from		"
" scale=.0001	scale factor to apply to header field values		"
"									"
" inv=0		weight by header value			 		"
" 		=1 weight by inverse of header value	 		"
"									"
" Notes:							 	"
" This code is initially written with offset weighting in mind, but may	"
" be used for other, user-specified schemes.				"
"									"
" The rationale for this program is to correct for unwanted linear	"
" amplitude trends with offset prior to either CMP stacking or AVO work."
" The code has to be edited should other functions of a keyword be required."
"									"
" The default form of the weighting is to multiply the amplitudes of the"
" traces by a factor of:    ( a + b*keyword).				"
"									"
" If key2=  header field is  set then this program uses the weighting	"
" values read from that header field, instead. Note, that because most	"
" header fields are integers, the scale=.0001 permits 10001 in the header"
" to represent 1.0001.							"
"									"
" To see the list of available keywords, type:    sukeyword  -o  <CR>	"
"									"
;

/*  Sun Feb 24 13:30:07 2013
  Automatically modified for usage in SeaSeis  */
namespace suweight {


/* Credits:
 * Author: CWP: John Stockwell  February 1999.
 * Written for Chris Walker of UniqueStep Ltd., Bedford, U.K.
 * inv option added by Garry Perratt (Geocon).
 *
 * header fields accessed: ns, keyword
 */

/**************** end self doc ***********************************************/


/* Globals */
segy tr;

void* main_suweight( void* args )
{

	char *key;	/* header key word from segy.h		*/
	char *type;     /* ... its type				*/
	int index;	/* ... its index			*/
	Value val;	/* ... its value			*/
	float fval;     /* ... its value cast to float		*/

	char *key2=NULL;      /* header key word from segy.h	*/
	char *type2=NULL;     /* ... its type			*/
	int index2=0;	/* ... its index			*/
	Value val2;	/* ... its value			*/
	cwp_Bool is_key2=cwp_false;	/* is key2 set?		*/
	float scale=0.0001;	/* ... scale parameter		*/

	float a;	/* ... constant weighting parameter     */
	float b;	/* ... variable weighting parameter     */
	int ns;		/* number of samples per trace	  	*/

	int inv=0;	/* inverse weighting flag		*/

	/* Initialize */
	cseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;
	cseis_su::csSUTraceManager* cs2su = suArgs->cs2su;
	cseis_su::csSUTraceManager* su2cs = suArgs->su2cs;
	int argc = suArgs->argc;
	char **argv = suArgs->argv;
	cseis_su::csSUGetPars parObj;

	void* retPtr = NULL;  /*  Dummy pointer for return statement  */
	su2cs->setSUDoc( sdoc_suweight );
	if( su2cs->isDocRequestOnly() ) return retPtr;
	parObj.initargs(argc, argv);

	try {  /* Try-catch block encompassing the main function body */


	/* Get optional parameters */
	if (!parObj.getparstring("key", &key))	 key = "offset";
	if (!parObj.getparfloat("a",&a))		a = 1.0;
	if (!parObj.getparfloat("b",&b))		b = .0005;

	/* Get key type and index */
	type = hdtype(key);
	index = getindex(key);

	/* Get key2 type and index */
	if(parObj.getparstring("key2", &key2)) {
		is_key2=cwp_true;
		
		/* Get key2 type and index */
		type2 = hdtype(key2);
		index2 = getindex(key2);
		if (!parObj.getparfloat("scale",&scale))	scale = .0001;
	}

	/* Get inversion parameter */
	if (!parObj.getparint("inv",&inv))		inv = 0;

        parObj.checkpars();

	/* Get info from first trace */
	if (!cs2su->getTrace(&tr))  throw cseis_geolib::csException("can't get first trace");
	ns = tr.ns;

	/* Loop through traces */
	do {
		register int i = 0;	     /* counter */

		if (!is_key2) { 
			/* Get value of key and convert to float */
			gethval(&tr, index, &val);
			fval = vtof(type,val);

			if (inv == 0) {
				/* Loop over samples in trace and apply weighting */
				for (i=0; i < ns; ++i)
					tr.data[i] *= (a + fval*b);
			} else {
				for (i=0; i < ns; ++i)
					tr.data[i] /= (a + fval*b);
			}
		} else {
			/* Get the value of key2 and convert to float */
			gethval(&tr, index2, &val2);
			fval = vtof(type2,val2);
			fval *= scale;

			if (inv == 0) {
				/* Loop over samples in trace and apply weighting */
				for (i=0; i < ns; ++i)
					tr.data[i] *= fval;
			} else {
				/* Loop over samples in trace and apply weighting */
				for (i=0; i < ns; ++i)
					tr.data[i] /= fval;
			}
		}

		/* Put out weighted traces */
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
