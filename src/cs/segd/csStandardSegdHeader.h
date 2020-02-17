/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_STANDARD_SEGD_HEADER_H
#define CS_STANDARD_SEGD_HEADER_H

namespace cseis_geolib {
  template<typename T> class csVector;
  class csHeaderInfo;
}

namespace cseis_segd {

/**
* Standard Segd Headers
* @author Bjorn Olofsson
*/
class csStandardSegdHeader {
public:
  static int const HDR_MAPPING_OBC      = 1111;
  static int const HDR_MAPPING_TOWED    = 1112;
  static int const HDR_MAPPING_STANDARD = 1113;

  static void setStandardHeaders( int segdHeaderMap, cseis_geolib::csVector<cseis_geolib::csHeaderInfo const*>* stdHdrList );
  ~csStandardSegdHeader();
private:
  csStandardSegdHeader();
};  // END class

}  // END namespace

#endif

