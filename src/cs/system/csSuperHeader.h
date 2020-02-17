/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SUPER_HEADER_H
#define CS_SUPER_HEADER_H

#include <string>
#include "geolib_defines.h"

namespace cseis_geolib {
  template<typename T> class csVector;
}

namespace cseis_system {

/**
 * Super header for Cseis seismic traces
 *
 * Stores headers that do not change from trace to trace, but may change from module to module.
 * Each module in Cseis flow has its unique super header
 * Super header defines super header values OUTPUT to the next module, i.e. the facts AFTER this module's exec phase
 * Sample interval/number of samples may be different for traces when INPUT to this module, see csTrace->getNumSamples()
 * For example see module TIME_GATE
 *
 * - sample interval
 * - Number of samples
 * - Ensemble keys
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
class csSuperHeader {
private:
  struct keyInfoStruct {
    std::string           name;
    cseis_geolib::type_t  hdrType;
    int                   hdrIndex;
  };
public:
  csSuperHeader();
  ~csSuperHeader();
  void set( csSuperHeader const* superHeader );

// Ensemble keys:
  void clearEnsembleKeys();
  /**
  * Set ensemble key
  * @param key       Name of header
  * @param type      Header type (INT, FLOAT..)
  * @param hdrIndex  Index of header (where to find it in header definition)
  */
  void setEnsembleKey( std::string const& keyName, int hdrIndex, cseis_geolib::type_t type, int keyIndex );
  void setEnsembleKey( std::string const& keyName, int keyIndex = 0 );
  int numEnsembleKeys() const;
  std::string const* ensembleKey( int keyIndex ) const;
  int ensembleKeyHeaderIndex( int keyIndex ) const {
    return myEnsKeyInfo[keyIndex].hdrIndex;
  }
  cseis_geolib::type_t ensembleKeyHeaderType( int keyIndex ) const {
    return myEnsKeyInfo[keyIndex].hdrType;
  }
  void removeEnsembleKey( int keyIndex );

  int mpi_getByteSize() const;
  int mpi_compress( char* data ) const;
  void mpi_decompress( char const* data );

  void dump( FILE* fout ) const;
 public:
  /// Number of samples in trace, including sample at time 0
  int   numSamples;
  /// Sample interval. Unit: [ms] or [Hz] or [m]
  float sampleInt;
  /// Domain, e.g. X-T, F-X, F-K...
  int   domain;
  /// In case of frequency domain, this contains the original number of samples in the time domain 
  int   numSamplesXT;
  /// In case of frequency domain, this contains the original sample interval [ms] in the time domain
  float sampleIntXT;
  /// Flag defining data type in frequency domain, e.g. amp_phase, real_imag...
  int   fftDataType;
  /// Survey grid definition: Grid origin X coordinate
  double grid_orig_x;
  /// Survey grid definition: Grid origin Y coordinate
  double grid_orig_y;
  /// Survey grid definition: Grid origin inline number (default = 1)
  int grid_orig_il;
  /// Survey grid definition: Grid origin xline number (default = 1)
  int grid_orig_xl;
  /// Survey grid definition: Bin size in inline direction (row/col step = 1)
  double grid_binsize_il;
  /// Survey grid definition: Bin size in xline direction (row/col step = 1)
  double grid_binsize_xl;
  /// Survey grid definition: Direction of inlines, clock-wise from North.
  double grid_azim_il;
  /// Survey grid definition: Direction of xlines, default is at -90deg from inline directon (left-hand coordinate system)
  double grid_azim_xl;
private:
  keyInfoStruct* myEnsKeyInfo;
  int myNumEnsKeys;
  int myNumAllocatedEnsKeys;
  void reallocate( int newNumAllocatedKeys );
  csSuperHeader( csSuperHeader const& obj );
};

} // namespace

#endif

