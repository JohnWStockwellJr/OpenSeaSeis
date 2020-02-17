/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEISMIC_READER_VER02_H
#define CS_SEISMIC_READER_VER02_H

#include <cstdio>
#include <string>
#include <fstream>
#include "geolib_defines.h"
#include "csSeismicReader_ver.h"


namespace cseis_io {

/**
 * Seismic file Reader, Cseis format
 *
 * Version 0.2
 * 
 * @author Bjorn Olofsson
 * @date 2010
 */
class csSeismicReader_ver02 : public csSeismicReader_ver {
 public:
  static int const VERSION_SEISMIC_READER   = 02;

 public:
  csSeismicReader_ver02( std::string filename, bool enableRandomAccess, int numTracesBuffer = 0 );
  virtual ~csSeismicReader_ver02();
  virtual bool readFileHeader( csSeismicIOConfig* config );
};

} // end namespace
#endif

