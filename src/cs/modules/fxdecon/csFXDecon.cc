/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "csFXDecon.h"
#include "csException.h"
#include "csFFT.h"
#include "geolib_defines.h"

using namespace mod_fxdecon;
using namespace std;

csFXDecon::csFXDecon() {
  myFFT = NULL;
  myNumTraces = 0;

  fdata = NULL;
  fdataw = NULL;
  ffreq = NULL;
  rmatrix = NULL;

  fvector = NULL;
  sfreq = NULL;
  autocorr = NULL;
  rautoc = NULL;
  iautoc = NULL;
  gvector = NULL;

  info = NULL;
  ipvt = NULL;
  bufferRealImag = NULL;

 sfreqout = NULL;
 tidataw  = NULL;
 ttodataw = NULL;
}
csFXDecon::~csFXDecon() {
  if( myFFT != NULL ) {
    delete myFFT;
    myFFT = NULL;
  }
  freeMem( myNumTraces );
}
//--------------------------------------------------------------------------------
//
void csFXDecon::initialize( float sampleInt_ms, int numSamplesIn, mod_fxdecon::Attr const& attr ) {
  myFFT = new cseis_geolib::csFFT( numSamplesIn );
  myNumSamplesFFT = myFFT->numFFTValues();
  // !CHANGE!    Try to reduce the number of samples in FFT: Feed fftTool maximum length of window (numSamplesWinF??) instead of full trace length (numSamples)
  myNumFreq       = myFFT->numFreqValues(); // myNumSamplesFFT/2 + 1;
  myFreqStep_hz   = 1000.0/(float)(myNumSamplesFFT*sampleInt_ms);
  mySampleInt_s   = sampleInt_ms / 1000.0;

  fmin           = attr.fmin;
  fmax           = attr.fmax;
  winLen_samp    = attr.winLen_samp;
  taperLen_samp  = attr.taperLen_samp;
  ntraces_design = attr.ntraces_design;
  ntraces_filter = attr.ntraces_filter;
  numWin         = attr.numWin;
  taperLen_s     = attr.taperLen_s;

  myNumSamplesWinF = winLen_samp + taperLen_samp/2;
  myNumSamplesWinI = winLen_samp + taperLen_samp;

  
  ntraces_fdataw = 2*ntraces_design+2*ntraces_filter;
}
void csFXDecon::freeMem( int numTraces ) {
  for( int itrc = 0; itrc < numTraces; itrc++ ) {
    delete [] fdata[itrc];
  }

  delete [] fdata;
  delete [] sfreq;
  delete [] autocorr;
  delete [] rautoc;
  delete [] gvector;
  delete [] fvector;

  for( int itrc = 0; itrc < 2*ntraces_filter; itrc++ ) {
    delete [] rmatrix[itrc];
  }
  delete [] rmatrix;
  delete [] info;
  delete [] ipvt;
  delete [] sfreqout;
  for( int itrc = 0; itrc < 2*ntraces_design; itrc++ ) {
    delete [] ffreq[itrc];
  }
  delete [] ffreq;
  delete [] ttodataw;
  delete [] bufferRealImag;
  delete [] iautoc;
  delete [] tidataw;
  for( int itrc = 0; itrc < ntraces_fdataw; itrc++ ) {
    delete [] fdataw[itrc];
  }
  delete [] fdataw;

  fdata = NULL;
  fdataw = NULL;
  ffreq = NULL;
  rmatrix = NULL;

  fvector = NULL;
  sfreq = NULL;
  autocorr = NULL;
  rautoc = NULL;
  iautoc = NULL;
  gvector = NULL;

  info = NULL;
  ipvt = NULL;
  bufferRealImag = NULL;

  sfreqout = NULL;
  tidataw = NULL;
  ttodataw = NULL;
}
//--------------------------------------------------------------------------------
//
void csFXDecon::initialize_internal( int numTraces, int numSamples ) {
  if( fdata != NULL ) {
    if( myNumTraces == numTraces ) {
      return;
    }
    freeMem( myNumTraces );
  }
  myNumTraces = numTraces;

  fdata  = new complex<float>*[ myNumTraces ];
  fdataw = new complex<float>*[ ntraces_fdataw ];
  ffreq  = new complex<float>*[ 2*ntraces_design ];
  rmatrix         = new float*[ 2*ntraces_filter ];

  for( int itrc = 0; itrc < numTraces; itrc++ ) {
    fdata[itrc]      = new complex<float>[myNumFreq];
  }
  for( int itrc = 0; itrc < ntraces_fdataw; itrc++ ) {
    fdataw[itrc]     = new complex<float>[myNumFreq];
  }
  for( int itrc = 0; itrc < 2*ntraces_filter; itrc++ ) {
    rmatrix[itrc]    = new float[ 2*ntraces_filter ];
  }
  for( int itrc = 0; itrc < 2*ntraces_design; itrc++ ) {
    ffreq[itrc] = new std::complex<float>[ myNumFreq ];
  }

  fvector  = new complex<float>[ 2*ntraces_filter+1 ];
  sfreq   = new complex<float>[2*ntraces_design+2*ntraces_filter];
  for( int itrc = 0; itrc < 2*ntraces_design+2*ntraces_filter; itrc++ ) {
    sfreq[itrc] = std::complex<float>(0,0);
  }
  autocorr = new complex<float>[ntraces_filter+1];
  rautoc  = new float[ ntraces_filter+1 ];
  iautoc  = new float[ ntraces_filter+1 ];
  gvector = new float[ 2*ntraces_filter ];

  info = new float[ 4*ntraces_design ];
  ipvt  = new int[ 4*ntraces_design ];
  bufferRealImag = new float[ 2*myNumSamplesFFT ];

  sfreqout = new std::complex<float>[ 2*ntraces_design ];
  tidataw  = new float[ numSamples ];
  ttodataw = new float[ myNumSamplesFFT ];
}
//--------------------------------------------------------------------------------
//
void csFXDecon::apply( float** samplesIn, float** samplesOut, int numTraces, int numSamples ) {
  initialize_internal( numTraces, numSamples );

  int numWinSpatial = numTraces / ntraces_design;

  int numSamplesWinCurrent = 0;

  // Loop over time windows
  for( int iwin = 0; iwin < numWin; iwin++ ) {

    // Zero arrays:
    for( int itrc = 0; itrc < numTraces; itrc++ ) {
      for( int ifq = 0; ifq < myNumFreq; ifq++ ) {
        fdata[itrc][ifq] = std::complex<float>(0,0);
      }
    }
    for( int itrc = 0; itrc < ntraces_design+2*ntraces_filter; itrc++ ) {
      for( int ifq = 0; ifq < myNumFreq; ifq++ ) {
        fdataw[itrc][ifq] = std::complex<float>(0,0);
      }
    }

    if( iwin > 0 && iwin < numWin-1 ) {
      numSamplesWinCurrent = myNumSamplesWinI;
    }
    else if ( iwin == 0 ) {
      if( numWin > 1 ) {
        numSamplesWinCurrent = myNumSamplesWinF;
      }
      else {
        numSamplesWinCurrent = numSamples;
      }
    }
    else {
      numSamplesWinCurrent = numSamples - winLen_samp*iwin + taperLen_samp/2;
    }

    // select data
    for( int itrc = 0; itrc < numTraces; itrc++ ) {
      if( iwin > 0 ) {
        memcpy( tidataw, &samplesIn[itrc][iwin*winLen_samp - taperLen_samp/2], numSamplesWinCurrent*sizeof(float) );
      }
      else {
        memcpy( tidataw, &samplesIn[itrc][0], numSamplesWinCurrent*sizeof(float) );
      }
      if( numSamples > numSamplesWinCurrent ) {
        memset( &tidataw[numSamplesWinCurrent], 0, (numSamples-numSamplesWinCurrent)*sizeof(float) );
      }

      myFFT->forwardTransform( tidataw, bufferRealImag, cseis_geolib::FX_REAL_IMAG );
      float const* realPtr = &bufferRealImag[0];
      float const* imagPtr = &bufferRealImag[myNumFreq];
      for( int ifq = 0; ifq < myNumFreq; ifq++ ) {
        fdata[itrc][ifq] = std::complex<float>( realPtr[ifq], imagPtr[ifq] );
        // fdata[itrc][ifq] = std::complex<float>( (float)realPtr[ifq] ,(float)imagPtr[ifq] );
        //  float amp = sqrt( fdata[itrc][ifq].real()*fdata[itrc][ifq].real() + fdata[itrc][ifq].imag()*fdata[itrc][ifq].imag() );
        //  float phase = atan2( -fdata[itrc][ifq].imag(), fdata[itrc][ifq].real() );
        //  fprintf(stdout,"%d %f  %f  %f  %f  %f\n",itrc,ifq*myFreqStep_hz,realPtr[ifq],imagPtr[ifq],amp,phase);
      }
    } // END for itrc

    //----------------------------------------------------------------------
    // Loop over space windows

    int ntrwu = 0;
    for( int jx = 0; jx < numWinSpatial; jx++ ) {

      // to take care of a possible incomplete last window
      if( numTraces < jx*ntraces_design+2*ntraces_design ) {
        ntrwu = numTraces - jx*ntraces_design;
      }
      else {
        ntrwu = ntraces_design;
      }

      // Select data
      for( int itrc = 0; itrc < ntrwu+2*ntraces_filter; itrc++ ) {
        for( int ifq = 0; ifq < myNumFreq; ifq++) {
          if( jx > 0 && jx < numWinSpatial-1 ) {
            fdataw[itrc][ifq] = fdata[itrc + jx*ntraces_design - ntraces_filter][ifq];
          }
          else if (jx==0) {
            if( itrc >= ntraces_filter && itrc < ntraces_design+ntraces_filter ) {
              fdataw[itrc][ifq] = fdata[itrc - ntraces_filter][ifq];
            }
            else if( itrc < ntraces_filter ) {
              fdataw[itrc][ifq] = fdata[0][ifq];
            }
            else {
              if( numWinSpatial > 1 ) {
                fdataw[itrc][ifq] = fdata[itrc - ntraces_filter][ifq];
              }
              else {
                fdataw[itrc][ifq] = fdata[numTraces-1][ifq];
              }
            }
          }
          else {
            if( itrc < ntrwu+ntraces_filter ) {
              fdataw[itrc][ifq] = fdata[itrc + jx*ntraces_design - ntraces_filter][ifq];
            }
            else {
              fdataw[itrc][ifq] = fdata[numTraces-1][ifq];
            }
          }
        } // for ifq
      } // for itrc


      // loop over frequencies
      for( int ifq = 0; ifq < myNumFreq; ifq++ ) {
        float freqCurrent = (float)ifq*myFreqStep_hz;
        if( freqCurrent >= fmin && freqCurrent <= fmax ) {
          
          // sfreq   = new complex<float>[2*ntraces_design+2*ntraces_filter];

          // Loop over space window
          for( int itrc = 0; itrc < ntrwu+2*ntraces_filter; itrc++ ) {
            sfreq[itrc]   = fdataw[itrc][ifq];
            //            fprintf(stdout,"%d %d  %f %d %f\n", iwin, jx, ifq*myFreqStep_hz, itrc, sfreq[itrc].real() );
          }
          
          // complex autocorrelation
          cxcor( ntrwu, 0, sfreq,   ntrwu, 0, sfreq,   ntraces_filter+1, 0, autocorr );
          
          // taking real and imaginary parts
          for( int itrc = 0; itrc < ntraces_filter+1; itrc++ ) {
            rautoc[itrc] = autocorr[itrc].real();
            iautoc[itrc] = autocorr[itrc].imag();
            //            fprintf(stdout,"%f %d %f %f\n", ifq*myFreqStep_hz, itrc, rautoc[itrc], iautoc[itrc]);
          }

          // zeroing files
          memset( (void*)gvector, 0, 2*ntraces_filter*sizeof(float) );
          for( int itrc = 0; itrc < 2*ntraces_filter+1; itrc++ ) {
            fvector[itrc] = std::complex<float>(0,0);
          }
          for( int itrc = 0; itrc < 2*ntraces_filter; itrc++ ) {
            memset( (void*)rmatrix[itrc], 0, 2*ntraces_filter*sizeof(float) );
          }

          // matrix problem
          for( int itrc = 0; itrc < ntraces_filter; itrc++ ) {
            for( int jtrc = 0; jtrc < ntraces_filter; jtrc++ ) { 
              if( itrc >= jtrc ) rmatrix[itrc][jtrc] = autocorr[itrc-jtrc].real();
              else        rmatrix[itrc][jtrc] = autocorr[jtrc-itrc].real();
            }
          }

          for( int itrc = ntraces_filter;itrc<2*ntraces_filter;itrc++) {
            for( int jtrc = 0;jtrc<ntraces_filter;jtrc++) {
              if( itrc-ntraces_filter < jtrc ) rmatrix[itrc][jtrc] = -autocorr[jtrc-itrc+ntraces_filter].imag();
              else            rmatrix[itrc][jtrc]= autocorr[itrc-jtrc-ntraces_filter].imag();
            }
          }

          for( int itrc = ntraces_filter; itrc < 2*ntraces_filter; itrc++ ) {
            for( int jtrc = ntraces_filter; jtrc < 2*ntraces_filter; jtrc++ )
              rmatrix[itrc][jtrc]=rmatrix[itrc-ntraces_filter][jtrc-ntraces_filter];
          }
          for( int itrc = 0;itrc<ntraces_filter;itrc++) {
            for( int jtrc = ntraces_filter; jtrc < 2*ntraces_filter; jtrc++ )
              rmatrix[itrc][jtrc] = -rmatrix[itrc+ntraces_filter][jtrc-ntraces_filter];
          }
          for( int itrc = 0; itrc < 2*ntraces_filter; itrc++ ) {
            if( itrc < ntraces_filter ) gvector[itrc] = autocorr[itrc+1].real();
            else gvector[itrc] = autocorr[itrc-ntraces_filter+1].imag();
          }

          csFXDecon::LUDecomposition( rmatrix, 2*ntraces_filter, ipvt, info );
          csFXDecon::LUBackSub( rmatrix, 2*ntraces_filter, ipvt, gvector );

          //for( int ii = 0; ii < 2*ntraces_filter; ii++ ) {
          //   for( int jj = 0; jj < 2*ntraces_filter; jj++ ) {
          //     fprintf(stdout,"%d %f %d %d %d %f  %f\n", jx, ifq*myFreqStep_hz, ii, jj, ipvt[ii*2*ntraces_filter+jj], info[ii*2*ntraces_filter+jj], rmatrix[ii][jj]);
          //   }
          //  }

          /* construct filter */
          for( int ifv = 0, ig = ntraces_filter-1; ifv < ntraces_filter; ifv++, ig-- ) {
            fvector[ifv] = std::conj( std::complex<float>( 0.5*gvector[ig], 0.5*gvector[ig+ntraces_filter] ) );
          }

          for( int ifv = ntraces_filter+1,ig=0 ;ifv < 2*ntraces_filter+1; ifv++,ig++ ) {
            fvector[ifv] = std::complex<float>( 0.5*gvector[ig], 0.5*gvector[ig+ntraces_filter] );
          }
          for( int itrc = 0; itrc < ntrwu; itrc++ ) {
            sfreqout[itrc] = std::complex<float>(0,0);
          }

          //for( int ifv = 0; ifv < 2*ntraces_filter+1; ifv++ ) {
          //  fprintf(stdout,"%d %f %d  %f %f\n", jx, ifq*myFreqStep_hz, ifv, fvector[ifv].real(), fvector[ifv].imag());
          //  }

          // convolution of data with filter
          // output is one sample ahead
          cconv( ntrwu+2*ntraces_filter, -ntraces_filter, sfreq,
                 2*ntraces_filter+1, -ntraces_filter, fvector,
                 ntrwu, 0, sfreqout ); 

          // store filtered values
          for( int itrc = 0; itrc < ntrwu; itrc++ ) {
            ffreq[itrc][ifq] = sfreqout[itrc];
            // fprintf(stdout,"%d %d %f  %f %f\n", jx, itrc, ifq*myFreqStep_hz, ffreq[itrc][ifq].real(), ffreq[itrc][ifq].imag());
          }

        }
      } // END for ifq frequencies loop


      // Loop along space windows
      for( int itrc = 0; itrc < ntrwu; itrc++ ) {
    
        for( int isamp = 0; isamp < myNumFreq; isamp++ ) {
          bufferRealImag[isamp]           = ffreq[itrc][isamp].real();
          bufferRealImag[isamp+myNumFreq] = ffreq[itrc][isamp].imag();
        }
        myFFT->inverseTransform( bufferRealImag, ttodataw, cseis_geolib::FX_REAL_IMAG );

        int minNumSamples = std::min(numSamples,myNumSamplesFFT);

        if( numSamples != myNumSamplesFFT ) {
          memset( &ttodataw[minNumSamples], 0, (std::max(numSamples,myNumSamplesFFT)-minNumSamples)*sizeof(float) );
        }

        // Loop along time
        if( numWin > 1 ) {
          // first portion of time window
          if( iwin > 0 ) {
            for( int isamp = 0; isamp < taperLen_samp; isamp++ ) {
              samplesOut[jx*ntraces_design+itrc][isamp+iwin*winLen_samp-taperLen_samp/2] +=
                ttodataw[isamp] * ( (float)isamp * mySampleInt_s / taperLen_s );
            }
          }
          else {
            for( int isamp = 0; isamp < taperLen_samp; isamp++ ) {
              samplesOut[jx*ntraces_design+itrc][isamp] = ttodataw[isamp];
            }
          }
          // intermediate portion of time window
          if( iwin > 0 ) {
            for( int isamp = taperLen_samp; isamp < numSamplesWinCurrent - taperLen_samp; isamp++ )
              samplesOut[jx*ntraces_design+itrc][isamp+iwin*winLen_samp-taperLen_samp/2] = ttodataw[isamp];
          }
          else {
            for( int isamp = taperLen_samp; isamp < numSamplesWinCurrent-taperLen_samp; isamp++ )
              samplesOut[jx*ntraces_design+itrc][isamp] = ttodataw[isamp];
          }
          // last portion of time window
          if( iwin > 0 && iwin < numWin-1 ) {
            for( int isamp = numSamplesWinCurrent-taperLen_samp; isamp < numSamplesWinCurrent; isamp++ )
              samplesOut[jx*ntraces_design+itrc][isamp+iwin*winLen_samp-taperLen_samp/2] +=
                ttodataw[isamp] * (1.0-((float)(isamp-numSamplesWinCurrent+taperLen_samp)) * mySampleInt_s / taperLen_s );
          }
          else if( iwin == numWin-1 ) {
            for( int isamp = numSamplesWinCurrent-taperLen_samp; isamp < numSamplesWinCurrent; isamp++ )
              samplesOut[jx*ntraces_design+itrc][isamp+iwin*winLen_samp-taperLen_samp/2] = ttodataw[isamp];
          }
          else {
            for( int isamp = numSamplesWinCurrent - taperLen_samp; isamp < numSamplesWinCurrent; isamp++ ) {
              samplesOut[jx*ntraces_design+itrc][isamp] += ttodataw[isamp] * ( 1.0 - ( (float)(isamp-numSamplesWinCurrent+taperLen_samp) ) * mySampleInt_s / taperLen_s );
            }
          }
        } // END if numWin > 1
        else {
          for( int isamp = 0; isamp < numSamples; isamp++ ) {
            samplesOut[jx*ntraces_design+itrc][isamp]=ttodataw[isamp];
          }
        }

      } // END for itrc - space windows

    } // END for jx - space windows
  } // END for iwin

}
void csFXDecon::dump() const {
  fprintf(stderr," --- FXDecon DUMP ---\n");
  fprintf(stderr,"numSamplesFFT:  %d\n", myNumSamplesFFT ); // nfft
  fprintf(stderr,"numFreq:        %d\n", myNumFreq );       // nf
  fprintf(stderr,"freqStep_hz:    %f\n", myFreqStep_hz );   // d1
  fprintf(stderr,"fmin:           %f\n", fmin );            // fmin
  fprintf(stderr,"fmax:           %f\n", fmax );            // fmax
  fprintf(stderr,"ntraces_design: %d\n", ntraces_design );  // ntrw
  fprintf(stderr,"ntraces_filter: %d\n", ntraces_filter );  // ntri / ntrf
  fprintf(stderr,"numWin:         %d\n", numWin );          // ntw
  fprintf(stderr,"winLen_samp:    %d\n", winLen_samp );     // nspws (twlen)
  fprintf(stderr,"taperLen_samp:  %d\n", taperLen_samp);    // nstaper
  fprintf(stderr,"taperLen_s:     %f\n", taperLen_s);       // taper
  fprintf(stderr,"numSamplesWinF: %d\n", myNumSamplesWinF );// nspwf
  fprintf(stderr,"numSamplesWinI: %d\n", myNumSamplesWinI );// nspwi
}

// complex correlation
void csFXDecon::cxcor( int num1, int index1, std::complex<float>* in1,
                       int num2, int index2, std::complex<float>* in2, 
                       int numCorr, int indexCorr, std::complex<float>* corr )
{
  complex<float>* xr = new complex<float>[ num1 ];
  for( int i = 0, j = num1-1; i < num1; ++i,--j ) {
    xr[i] = std::conj( in1[j] );
  }
  csFXDecon::cconv( num1, 1-index1-num1, xr, num2, index2, in2,  numCorr, indexCorr, corr );
  delete [] xr;
}

// complex convolution
void csFXDecon::cconv( int num1, int index1, std::complex<float>* in1,
                       int num2, int index2, std::complex<float>* in2, 
                       int numCorr, int indexCorr, std::complex<float>* corr )
{
  int ilx = index1+num1-1;
  int ily = index2+num2-1;
  int ilz = indexCorr+numCorr-1;

  for( int i = indexCorr; i <= ilz; ++i ) {
    int jlow = i-ily;
    if( jlow < index1 ) jlow = index1;
    int jhigh = i-index2;
    if( jhigh > ilx ) jhigh = ilx;
    complex<float> sum = std::complex<float>( 0, 0) ;
    for( int j = jlow; j <= jhigh; ++j ) {
      sum += in1[j-index1] * in2[i-j-index2];
    }
    corr[i-indexCorr] = sum;
  }
}


void csFXDecon::LUDecomposition( float** AA, int nIn, int* indx, float* dd )
{
  float dum, sum, temp;
  unsigned int imax = 0;
  unsigned int nIn_unsigned = (unsigned int)nIn;

  float* vv = new float[ nIn_unsigned ];
  if( !vv ) fprintf(stderr,"Allocation failure in ludcmpBO\n");

  *dd = 1.0;
  for( unsigned int i = 0; i < nIn_unsigned; i++ ) {
    float big = 0.0;
    for( unsigned int j = 0; j < nIn_unsigned; j++ ) {
      if ((temp = fabs(AA[i][j])) > big) big = temp;
    }
    if (big == 0.0) fprintf(stderr,"Singular matrix in routine LUDecomposition\n");
    vv[i] = 1.0 / big;
  }

  for( unsigned int j = 0; j < nIn_unsigned; j++ ) {
    for( unsigned int i = 0; i < j; i++ ) {
      sum = AA[i][j];
      for( unsigned int k = 0; k < i; k++) sum -= AA[i][k]*AA[k][j];
      AA[i][j] = sum;
    }
    float big = 0.0;
    for( unsigned int i = j; i < nIn_unsigned; i++ ) {
      sum = AA[i][j];
      for( unsigned int k = 0; k < j; k++ ) {
        sum -= AA[i][k]*AA[k][j];
      }
      AA[i][j] = sum;
      if( (dum = vv[i]*fabs(sum)) >= big ) {
        big = dum;
        imax = i;
      }
    }
    if( j != imax ) {
      for( unsigned int k = 0; k < nIn_unsigned; k++ ) {
        dum = AA[imax][k];
        AA[imax][k] = AA[j][k];
        AA[j][k] = dum;
      }
      *dd = -(*dd);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if( AA[j][j] == 0.0 ) AA[j][j] = 1.0e-10;
    if( j != nIn_unsigned ) {
      dum = 1.0 / AA[j][j];
      for( unsigned int i = j+1; i < nIn_unsigned; i++ ) AA[i][j] *= dum;
    }
  }

  delete [] vv;
}

void csFXDecon::LUBackSub( float** AA, int nIn, int* indx, float* bb )
{
  int ii = -1;

  for( int i = 0; i < nIn; i++ ) {
    int ip = indx[i];
    float sum = bb[ip];
    bb[ip] = bb[i];
    if( ii != -1 ) {
      for( int j = ii; j <= i-1; j++ ) sum -= AA[i][j]*bb[j];
    }
    else if( sum ) {
      ii = i;
    }
    bb[i] = sum;
  }
  for( int i = nIn-1; i >= 0; i-- ) {
    float sum = bb[i];
    for( int j = i+1; j < nIn; j++ ) sum -= AA[i][j]*bb[j];
    bb[i] = sum/AA[i][i];
  }
}
