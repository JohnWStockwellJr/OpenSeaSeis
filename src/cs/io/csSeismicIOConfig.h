/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEISMIC_IO_CONFIG_H
#define CS_SEISMIC_IO_CONFIG_H

#include "csVector.h"
#include "geolib_defines.h"
#include <string>

namespace cseis_geolib {
  class csHeaderInfo;
}

namespace cseis_io {

/**
 * Configuration for Cseis seismic reader/writer
 * Contains all necessary parameters information to facilitate reading/writing of seismic files, i.e.
 * super header information such as nmber of samples, and trace header information
 *
 * @author Bjorn Olofsson
 * @date 2007
 */
class csSeismicIOConfig {
 public:
  csSeismicIOConfig();
  ~csSeismicIOConfig();
  /**
  * Add information about trace header
  */
  void addHeader( cseis_geolib::csHeaderInfo const* info );
  /**
  * Add information about trace header
  */
  void addHeader( cseis_geolib::type_t type, std::string const& name, std::string const& description, int nElements );
  /**
  * @return header info for given trace header
  */
  cseis_geolib::csHeaderInfo const* headerInfo( int headerIndex ) const;
  cseis_geolib::csHeaderInfo const* headerInfo( std::string const& headerName ) const;
  int headerIndex( std::string const& headerName ) const;  
  
  /**
  * @return number of trace headers
  */
  int numTrcHeaders() const;
  /**
  * @param keyIndex (i) Index of ensemble key
  * @return name of given ensemble key
  */
  std::string const* keyName( int keyIndex );
  /**
  *
  */
  int numEnsKeys() const;
    
public:
  // Super header
  int     numSamples;
  float   sampleInt;
  int     domain;
  int     fftDataType;
  int     numSamplesXT;
  float   sampleIntXT;
  double  grid_orig_x;
  double  grid_orig_y;
  int     grid_orig_il;
  int     grid_orig_xl;
  double  grid_binsize_il;
  double  grid_binsize_xl;
  double  grid_azim_il;
  double  grid_azim_xl;
  cseis_geolib::csVector<std::string> ensKeyNames;
  int     numTraces;

  // Trace header
  int     byteSizeSamples;  // Number of bytes in one trace (samples only)
  int     byteSizeHdrValueBlock;  // Number of bytes in one trace header
private:
  cseis_geolib::csVector<cseis_geolib::csHeaderInfo const*> myTrcHdrList;
};

} // namespace
#endif

