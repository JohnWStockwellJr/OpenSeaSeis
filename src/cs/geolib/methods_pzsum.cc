/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "geolib_math.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>

using namespace std;

namespace cseis_geolib {

bool Powerof2( int nx, int* m, int* twopm);
bool fft(int dir,int m,double *x,double *y);
void work_fft( float const* data_xt, double* bufferReal, double* bufferImag, int nSamplesXT, int nSamplesFFT, int twoPower, double* amp_spectrum  );
void work_fft2( float const* data_xt, double* bufferReal, double* bufferImag, int nSamplesXT, int nSamplesFFT, int twoPower );

/**
* PZ analysis in XT domain
*/
void pz_analysis( float const* const* samples_p,  // (i) Hydrophone samples, for all input traces
                  float const* const* samples_z,  // (i) Geophone samples, for all input traces
                  int    nTraces,           // (i) Number of input traces
                  int    nSamples,          // (i) Number of input samples
                  float* refl_coef,         // (o) Seabed reflection coefficient
                  float* pz_scalar,         // (o) PZ scalar
                  float* wdep,              // (o) Water depth [m]
                  float* z_phase )          // (o) Relative phase of geophone [rad]
{


  // (1) Compute auto-correlations in x-t  (should really be done in frequency domain...)
  int nSamplesCorr  = compute_correlation_length( nSamples/2 );
  float* p_auto_xt = new float[nSamplesCorr];
  float* z_auto_xt = new float[nSamplesCorr];
  float* pz_cross_xt = new float[nSamplesCorr];
  float* zp_cross_xt = new float[nSamplesCorr];
  float pzscal;
  
  if( !p_auto_xt || !z_auto_xt || !pz_cross_xt || !zp_cross_xt ) {
    fprintf(stderr,"Error: Out of memory...\n");
    return;
  }

  int maxlag = nSamples / 2;
  
  for( int itrc = 0; itrc < nTraces; itrc++ ) {
    compute_twosided_correlation( samples_p[itrc], samples_p[itrc], nSamples, p_auto_xt, maxlag );
    compute_twosided_correlation( samples_z[itrc], samples_z[itrc], nSamples, z_auto_xt, maxlag );
    compute_twosided_correlation( samples_p[itrc], samples_z[itrc], nSamples, pz_cross_xt, maxlag );
    compute_twosided_correlation( samples_z[itrc], samples_p[itrc], nSamples, zp_cross_xt, maxlag );

    pzscal = sqrt(p_auto_xt[nSamplesCorr/2] / z_auto_xt[nSamplesCorr/2]);
    fprintf(stderr,"%f %f  %f\n", p_auto_xt[nSamplesCorr/2], z_auto_xt[nSamplesCorr/2], pzscal);

    for( int i = 0; i < nSamplesCorr; i++ ) {
      //  fprintf(stdout,"%d  %f %f %f %f\n", i+1, p_auto_xt[i], pzscal*pzscal*z_auto_xt[i], pzscal*pz_cross_xt[i], pzscal*zp_cross_xt[nSamplesCorr-i-1] );
      //    fprintf(stdout,"%d  %f %f %f %f\n", i+1, p_auto_xt[i], pzscal*pzscal*z_auto_xt[i], pzscal*pz_cross_xt[i], pzscal*zp_cross_xt[i] );
    }
  }

  // (2) Compute fft's of all 4 auto-correlations
  int nSamplesFFT = nSamplesCorr;
  int twoPower;        // 2^twoPower: Number of samples for FFT
  int two_power_m;
  Powerof2( nSamplesCorr, &twoPower, &two_power_m );
  if( two_power_m != nSamplesCorr ) {
    nSamplesFFT  = two_power_m * 2;
    twoPower    += 1;
  }
  float sampleInt = 2.0;
  float df = 1000.0/((float)nSamplesFFT*sampleInt);

  double* p_fft_power = new double[nSamplesFFT];
  double* z_fft_power = new double[nSamplesFFT];
  double* pz_fft_power = new double[nSamplesFFT];
  double* zp_fft_power = new double[nSamplesFFT];
  
  double* bufferReal = new double[nSamplesFFT];
  double* bufferImag = new double[nSamplesFFT];

  work_fft( p_auto_xt, bufferReal, bufferImag, nSamplesCorr, nSamplesFFT, twoPower, p_fft_power  );
  work_fft( z_auto_xt, bufferReal, bufferImag, nSamplesCorr, nSamplesFFT, twoPower, z_fft_power  );
  work_fft( pz_cross_xt, bufferReal, bufferImag, nSamplesCorr, nSamplesFFT, twoPower, pz_fft_power  );
  work_fft( zp_cross_xt, bufferReal, bufferImag, nSamplesCorr, nSamplesFFT, twoPower, zp_fft_power  );
  
  //
  double noise = 0.1;  
   for( int i = 0; i < nSamplesFFT; i++ ) {
     fprintf(stdout,"%d  %f %f %f %f   %f  %f\n", i+1, p_fft_power[i], z_fft_power[i], pz_fft_power[i], zp_fft_power[i], df*i,
             (pz_fft_power[i] + noise) * (zp_fft_power[i] + noise) /  (p_fft_power[i]+noise) );
   }

  
  delete [] p_fft_power;
  delete [] z_fft_power;
  delete [] pz_fft_power;
  delete [] zp_fft_power;
  
  
  delete [] p_auto_xt;
  delete [] z_auto_xt;
  delete [] pz_cross_xt;
  delete [] zp_cross_xt;
}

  //}

//----------------------------------------------------------------------------
/**
* PZ analysis in XT domain
*/
void pz_analysis2( float const* const* p_samples,  // (i) Hydrophone samples, for all input traces
                   float const* const* z_samples,  // (i) Geophone samples, for all input traces
                   int    nTraces,           // (i) Number of input traces
                   int    nSamples,          // (i) Number of input samples
                   float  waterVel,         // (i) Water velocity [m/s]
                   float* refl_coef,         // (o) Seabed reflection coefficient
                   float* pz_scalar,         // (o) PZ scalar
                   float* wdep,              // (i/o) Water depth [m]
                   float* z_phase )          // (o) Relative phase of geophone [rad]
{
  // (1) Compute FFT of P and Z
  int nSamplesFFT = nSamples;
  int twoPower;        // 2^twoPower: Number of samples for FFT
  int two_power_m;
  Powerof2( nSamplesFFT, &twoPower, &two_power_m );
  if( two_power_m != nSamplesFFT ) {
    nSamplesFFT  = two_power_m * 2;
    twoPower    += 1;
  }
  float sampleInt = 2.0;
  float df = 1000.0/((float)nSamplesFFT*sampleInt);

  double* p_fft_real = new double[nSamplesFFT];
  double* z_fft_real = new double[nSamplesFFT];
  double* p_fft_imag = new double[nSamplesFFT];
  double* z_fft_imag = new double[nSamplesFFT];

  if( !p_fft_real || !p_fft_imag || !z_fft_real || !z_fft_imag ) {
    fprintf(stderr,"Error: Out of memory...\n");
    return;
  }

  fprintf(stderr,"nSamplesFFT: %d, Water depth: %f, Water vel: %f\n", nSamplesFFT, *wdep, waterVel);
  
  for( int itrc = 0; itrc < nTraces; itrc++ ) {
    work_fft2( p_samples[itrc], p_fft_real, p_fft_imag, nSamples, nSamplesFFT, twoPower );
    work_fft2( z_samples[itrc], z_fft_real, z_fft_imag, nSamples, nSamplesFFT, twoPower );

    
    //  for( int i = 0; i < nSamplesFFT; i++ ) {
    //   double freq = df * i;
    //  if( i > nSamplesFFT/2 ) freq = (i-nSamplesFFT)*df;
    //  fprintf(stdout,"%d %f   %f %f %f   %f %f %f\n", i+1, freq,
    //          p_fft_real[i], p_fft_imag[i], sqrt(pow(p_fft_real[i],2) + pow(p_fft_imag[i],2)),
    //          z_fft_real[i], z_fft_imag[i], sqrt(pow(z_fft_real[i],2) + pow(z_fft_imag[i],2)) );
    // }
    
    
    // (2) Scan for PZ scalar, minimisation
    complex<double> z;
    //    double omega;
    double pzscalMin = 5.0;
    double pzscalMax = 20.0;
    double minFreq   = 0.0;
    double maxFreq   = 70.0;
    double pzscalRange = pzscalMax - pzscalMin;
    int nSteps = 101;
    int nFreqs = nSamplesFFT;
    complex<double> z_value;
    complex<double> p_value;
    double waterDepth = *wdep + 0.0;
    double inciAngle  = 0.0;
    complex<double> i_complex( 0.0, 1.0 );

    //for( int istep = 0; istep < nSteps; istep++ ) {
    int istep = 0;
      double pzscal = pzscalMin + (double)istep/(double)(nSteps-1) * pzscalRange;
      double sum = 0.0;
      int counter = 0;
      for( int ifreq = 0; ifreq < nFreqs; ifreq++ ) {
        double freq = ifreq * df;
        if( freq < minFreq || freq > maxFreq ) continue;
        if( ifreq > nSamplesFFT/2 ) freq = (ifreq-nSamplesFFT)*df;
        z = exp( i_complex * 2.0*M_PI*freq * (2.0*waterDepth/waterVel) * cos(inciAngle) );
        p_value = complex<double>( p_fft_real[ifreq], p_fft_imag[ifreq] );
        z_value = complex<double>( z_fft_real[ifreq], z_fft_imag[ifreq] );
        complex<double> result = pzscal*(1.0 - z)*z_value - (1.0 + z)*p_value*cos(inciAngle);
        //complex<double> result = pzscal*z_value - p_value*cos(inciAngle);
        //double zr = complex<double>( (1.0 - z)*z_value ).real();
        //double zi = complex<double>( (1.0 - z)*z_value ).imag();
        //double pr = complex<double>( (1.0 + z)*p_value ).real();
        //double pi = complex<double>( (1.0 + z)*p_value ).imag();
        sum += abs(result);
        counter += 1;
        //   fprintf(stdout,"%d  %f   %f %f %f %f    %f %f %f %f   %f %f %f\n", istep, freq,
        //        pr, 5*zr, p_fft_real[ifreq], 5*z_fft_real[ifreq],
        //        pi, 5*zi, p_fft_imag[ifreq], 5*z_fft_imag[ifreq],
        //        complex<double>(1.0-z).real(), complex<double>(1.0-z).imag(), abs(1.0-z) );
        //fprintf(stdout,"%d  %f  real: %f %f %f   imag: %f %f %f  %f\n", istep, freq, pr, 5*zr, pr-5*zr, pi, 2.0*M_PI*freq, (2.0*waterDepth/waterVel), 2.0*M_PI*freq * (2.0*waterDepth/waterVel) * cos(inciAngle)  );
      }
      //           fprintf(stdout,"PZSCAL  %d %f %f\n", istep, pzscal, sum/counter);
      // }

      /*
    // (3) Scan for seabed reflectivity, minimisation
    double reflMin = 0.0;
    double reflMax = 0.9;
    double reflRange = reflMax - reflMin;
    double pzscal_final = 5.0;
    for( int istep = 0; istep < nSteps; istep++ ) {
      double refl = reflMin + (double)istep/(double)(nSteps-1) * reflRange;
      double sum = 0.0;
      int counter = 0;
      for( int ifreq = 0; ifreq < nFreqs; ifreq++ ) {
        double freq = ifreq * df;
        if( freq < minFreq || freq > maxFreq ) continue;
        if( ifreq > nSamplesFFT/2 ) freq = (ifreq-nSamplesFFT)*df;
        p_value = complex<double>( p_fft_real[ifreq], p_fft_imag[ifreq] );
        z_value = complex<double>( z_fft_real[ifreq], z_fft_imag[ifreq] );
        complex<double> result = pzscal_final*(1.0 + refl)*z_value + (1.0 - refl)*p_value*cos(inciAngle);
        sum += abs(result);
        counter += 1;
        // fprintf(stdout,"%d  %f   %f %f %f %f    %f %f %f %f   %f %f %f\n", istep, freq,
        //         p_value.real(), p_value.imag(),
        //        pi, 5*zi, p_fft_imag[ifreq], 5*z_fft_imag[ifreq],
        //        complex<double>(1.0-z).real(), complex<double>(1.0-z).imag(), abs(1.0-z) );
        // fprintf(stdout,"%d  %f  real: %f %f %f   imag: %f %f %f  %f\n", istep, freq, pr, 5*zr, pr-5*zr, pi, 2.0*M_PI*freq, (2.0*waterDepth/waterVel), 2.0*M_PI*freq * (2.0*waterDepth/waterVel) * cos(inciAngle)  );
      }
      fprintf(stdout,"REFL  %d %f %f\n", istep, refl, sum/counter);
    }
      */
  }  // END: for itrc

  
  delete [] p_fft_real;
  delete [] z_fft_real;
  delete [] p_fft_imag;
  delete [] z_fft_imag;
}

void work_fft( float const* data_xt, double* bufferReal, double* bufferImag, int nSamplesXT, int nSamplesFFT, int twoPower, double* amp_spectrum  ) {
  int direction= 1;
  for( int i = 0; i < nSamplesXT; i++ ) {
    bufferReal[i] = data_xt[i];
    bufferImag[i] = 0;
  }
  for( int i = nSamplesXT; i < nSamplesFFT; i++ ) {
    bufferReal[i] = 0;
    bufferImag[i] = 0;
  }
  if( !fft( direction, twoPower, bufferReal, bufferImag ) ) {
    fprintf(stderr,"FFT transform failed for unknown reasons...\n");
    exit(-1);
  }
  for( int i = 0; i < nSamplesFFT; i++ ) {
    amp_spectrum[i] = sqrt( bufferReal[i]*bufferReal[i] + bufferImag[i]*bufferImag[i] );
    //    fft_power[i] = bufferReal[i];
    //    fft_power[i] = bufferImag[i];
  }

}

void work_fft2( float const* data_xt, double* bufferReal, double* bufferImag, int nSamplesXT, int nSamplesFFT, int twoPower ) {
  int direction = 1;
  for( int i = 0; i < nSamplesXT; i++ ) {
    bufferReal[i] = data_xt[i];
    bufferImag[i] = 0;
  }
  for( int i = nSamplesXT; i < nSamplesFFT; i++ ) {
    bufferReal[i] = 0;
    bufferImag[i] = 0;
  }
  if( !fft( direction, twoPower, bufferReal, bufferImag ) ) {
    fprintf(stderr,"FFT transform failed for unknown reasons...\n");
    exit(-1);
  }
}




bool Powerof2( int nx, int* m, int* twopm) {
  int value = nx;
  *m = 0;
  while( (value = (int)(value / 2)) > 0 ) {
    *m += 1;
  }
  *twopm = (int)pow( 2.0, *m );
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


}
