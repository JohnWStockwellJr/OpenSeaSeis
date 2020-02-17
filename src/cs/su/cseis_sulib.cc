/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSUGetPars.h"
#include "csException.h"
#include "par.h"
#include "su.h"
#include "segy.h"

void cseis_decodeReflectors( cseis_su::csSUGetPars* parObj, int *nrPtr,
                             float **aPtr, int **nxzPtr, float ***xPtr, float ***zPtr)
/*************************************************************************
decodeReflectors - parse reflectors parameter string
**************************************************************************
Output:
nrPtr		pointer to nr an int specifying number of reflectors
aPtr		pointer to a specifying reflector amplitudes
nxzPtr		pointer to nxz specifying number of (x,z) pairs defining the
		reflectors
xPtr		pointer to array[x][nr] of x values for the entire model
zPtr		array[z][nr] of z values for the entire model

***************************************************************************
Author: Dave Hale, Colorado School of Mines, 09/17/91
**************************************************************************/
{
	int nr,*nxz,ir;
	float *a,**x,**z;
	char t[1024],*s;

	/* count reflectors */
	nr = parObj->countparname("ref");
	if (nr==0) nr = 1;
	
	/* allocate space */
	a   = (float*)ealloc1(nr,sizeof(float));
	nxz = (int*)ealloc1(nr,sizeof(int));
	x   = (float**)ealloc1(nr,sizeof(float*));
	z   = (float**)ealloc1(nr,sizeof(float*));

	/* get reflectors */
	for (ir=0; ir<nr; ++ir) {
		if (!parObj->getnparstring(ir+1,"ref",&s)) s = "1:1,2;4,2";
		strcpy(t,s);
		if (!decodeReflector(t,&a[ir],&nxz[ir],&x[ir],&z[ir]))
			err("Reflector number %d specified "
				"incorrectly!\n",ir+1);
	}

	/* set output parameters before returning */
	*nrPtr = nr;
	*aPtr = a;
	*nxzPtr = nxz;
	*xPtr = x;
	*zPtr = z;
}



/* Value getpar -- omitted string, bit types for now */
void cseis_getParVal( cseis_su::csSUGetPars* parObj, cwp_String name, cwp_String type, int n, Value *valp)
{
        register int k;
	short *h;
	unsigned short *u;
	long *l;
	unsigned long *v;
	int *i;
	unsigned int *p;
	float *f;
	double *d;
	
	switch(*type) {
        case 'h':
		h = (short*) ealloc1(n, sizeof(short));
		parObj->getparshort(name, h);  
		for (k = 0; k < n; ++k) valp[k].h = h[k];
	break;
        case 'u':
		u = (unsigned short*) ealloc1(n, sizeof(unsigned short));
		parObj->getparushort(name, u);  
		for (k = 0; k < n; ++k) valp[k].u = u[k];
	break;
        case 'i':
		i = (int*) ealloc1(n, sizeof(int));
		parObj->getparint(name, i);  
		for (k = 0; k < n; ++k) valp[k].i = i[k];
	break;  
        case 'p':
		p = (unsigned int*) ealloc1(n, sizeof(unsigned int));
		parObj->getparuint(name, p);  
		for (k = 0; k < n; ++k) valp[k].p = p[k];
	break;
        case 'l':
		l = (long*) ealloc1(n, sizeof(long));
		parObj->getparlong(name, l);  
		for (k = 0; k < n; ++k) valp[k].l = l[k];
	break;
        case 'v':
		v = (unsigned long*) ealloc1(n, sizeof(unsigned long));
		parObj->getparulong(name, v);  
		for (k = 0; k < n; ++k) valp[k].v = v[k];
	break;
        case 'f':
		f = (float*) ealloc1(n, sizeof(float));
		parObj->getparfloat(name, f);  
		for (k = 0; k < n; ++k) valp[k].f = f[k];
	break;  
        case 'd':
		d = (double*) ealloc1(n, sizeof(double));
		parObj->getpardouble(name, d);  
		for (k = 0; k < n; ++k) valp[k].d = d[k];
	break;  
        default:
                throw cseis_geolib::csException("cseis_getparval: %d: mysterious type %s", __LINE__, type);
        }
}
