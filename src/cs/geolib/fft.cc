/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "geolib_math.h"
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;

/**
* FFT transforms
*
*/

bool Powerof2( int nx, int* m, int* twopm) {
  int value = nx;
  *m = 0;
  while( (value = (int)(value / 2)) > 0 ) {
    *m += 1;
  }
  float mFloat = (float)*m;
  *twopm = (int)pow( 2.0f, mFloat );
  return( *twopm == nx );
}
bool fft( int dir, int power_of_two, double *realValues, double *imagValues );
bool fft( int dir, int power_of_two, double *realValues, double *imagValues, bool forceNormalisation );

/*
   This computes an in-place complex-to-complex FFT 
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform 
*/
bool fft_1d(short int dir,long m,double *x,double *y)
{
   long n,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   n = 1;
   for (i=0;i<m;i++) 
      n *= 2;

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) 
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }
   
   return(true);
}


/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there are memory problems or
      the dimensions are not powers of 2
*/
//bool FFT2D(complex<double>** c,int nx,int ny,int dir)

bool fft_2d(double** c_real, double** c_imag, int nx, int ny, int dir )
{
   int mx,my,twopm;
   double *real,*imag;

   // Transform the rows
   int max_nxny = cseis_geolib::MAX(nx,ny);
   real = new double[max_nxny];
   imag = new double[max_nxny];
   if (real == NULL || imag == NULL)
      return(false);
   if (!Powerof2(nx,&mx,&twopm) || twopm != nx)
      return(false);
   if (!Powerof2(ny,&my,&twopm) || twopm != ny)
      return(false);

   for(int j=0;j<ny;j++) {
      for(int i=0;i<nx;i++) {
         real[i] = c_real[i][j];
         imag[i] = c_imag[i][j];
      }
      fft(dir,mx,real,imag);
      for( int i=0; i<nx; i++ ) {
         c_real[i][j] = real[i];
         c_imag[i][j] = imag[i];
      }
   }

   // Transform the columns
  for( int i=0; i<nx; i++) {
//    memcpy( real, &c_real[i][0], ny*sizeof(double) );
//    memcpy( imag, &c_imag[i][0], ny*sizeof(double) );
      for( int j = 0; j < ny; j++ ) {
         real[j] = c_real[i][j];
         imag[j] = c_imag[i][j];
      }
    fft(dir,my,real,imag);
//    memcpy( &c_real[i][0], real, ny*sizeof(double) );
//    memcpy( &c_imag[i][0], imag, ny*sizeof(double) );
      for (int j=0;j<ny;j++) {
         c_real[i][j] = real[j];
         c_imag[i][j] = imag[j];
      }
   }
   delete [] real;
   delete [] imag;

   return(true);
}

/*-------------------------------------------------------------------------
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform

     Formula: forward
                  N-1
                  ---
              1   \          - j k 2 pi n / N
      X(n) = ---   >   x(k) e                    = forward transform
              N   /                                n=0..N-1
                  ---
                  k=0

      Formula: reverse
                  N-1
                  ---
                  \          j k 2 pi n / N
      X(n) =       >   x(k) e                    = forward transform
                  /                                n=0..N-1
                  ---
                  k=0
*/
bool fft( int dir, int power_of_two, double *realValues, double *imagValues ) {
  return fft( dir, power_of_two, realValues, imagValues, true );
}
bool fft( int dir, int power_of_two, double *realValues, double *imagValues, bool forceNormalisation )
{
   long nn,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   nn = 1;
   for (i=0;i<power_of_two;i++)
      nn *= 2;

   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn-1;i++) {
      if (i < j) {
         tx = realValues[i];
         ty = imagValues[i];
         realValues[i] = realValues[j];
         imagValues[i] = imagValues[j];
         realValues[j] = tx;
         imagValues[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<power_of_two;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<nn;i+=l2) {
            i1 = i + l1;
            t1 = u1 * realValues[i1] - u2 * imagValues[i1];
            t2 = u1 * imagValues[i1] + u2 * realValues[i1];
            realValues[i1] = realValues[i] - t1;
            imagValues[i1] = imagValues[i] - t2;
            realValues[i] += t1;
            imagValues[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for inverse transform */
   if( dir == -1 || forceNormalisation ) {
      for (i=0;i<nn;i++) {
         realValues[i] /= (double)nn;
         imagValues[i] /= (double)nn;
      }
   }

   return(true);
}

bool fft_float(int dir,int power_of_two,float *realValues,float *imagValues)
{
   long nn,i,i1,j,k,i2,l,l1,l2;
   float c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   nn = 1;
   for (i=0;i<power_of_two;i++)
      nn *= 2;

   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn-1;i++) {
      if (i < j) {
         tx = realValues[i];
         ty = imagValues[i];
         realValues[i] = realValues[j];
         imagValues[i] = imagValues[j];
         realValues[j] = tx;
         imagValues[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<power_of_two;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<nn;i+=l2) {
            i1 = i + l1;
            t1 = u1 * realValues[i1] - u2 * imagValues[i1];
            t2 = u1 * imagValues[i1] + u2 * realValues[i1];
            realValues[i1] = realValues[i] - t1;
            imagValues[i1] = imagValues[i] - t2;
            realValues[i] += t1;
            imagValues[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<nn;i++) {
         realValues[i] /= (float)nn;
         imagValues[i] /= (float)nn;
      }
   }

   return(true);
}

