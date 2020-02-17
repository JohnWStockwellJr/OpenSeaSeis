/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGY_HDR_MAP_H
#define CS_SEGY_HDR_MAP_H

#include <string>
#include <cstdio>
#include "geolib_defines.h"

namespace cseis_geolib {

class csHeaderInfo;
class csSegyHeaderInfo;
class csFlexHeader;
template<typename T> class csVector;

/**
* SEG-Y header mapping definition
*
* Defines byte mapping of SEG-Y trace headers
*
* @author Bjorn Olofsson
* @date 2006
*/
class csSegyHdrMap {
public:  
  static int const NONE = -1;
  static int const SEGY_STANDARD = 1110;
  static int const SEGY_OBC      = 1111;
  static int const SEGY_SEND     = 1112;
  static int const SEGY_ARMSS    = 1113;
  static int const SEGY_PSEGY    = 1114;
  static int const SEGY_NODE_OLD = 1115;
  static int const SEGY_NODE     = 1116;
  static int const SEGY_SU       = 1117;  // Special data format (Seismic Unix). Use this format to set Seaseis or SU standard header names (default)
  static int const SEGY_SU_ONLY  = 1118;  // Special data format (Seismic Unix). Use this format to set SU standard header names only
  static int const SEGY_SU_BOTH  = 1119;  // Special data format (Seismic Unix). Use this format to set Seaseis and SU standard header names (duplication)
  static int const SEGY_DAS      = 1120;

  static type_t const SEGY_HDR_TYPE_6BYTE = 255;
public:
  /**
  * Constructor
  * @param defaultMap  Default SEG-Y trace header mapping definition
  * @param initScalar  true if coordinate & elevation scalars shall be initialized immediately
  *                    If set to false, scalars must be initialized manually later on, by calling method initScalars().
  *                    See initScalars() for more information
  */
  csSegyHdrMap( int defaultMap, bool initScalars = true );
  csSegyHdrMap( int defaultMap, bool initScalars, std::string filename_hdrmap );
  csSegyHdrMap( csSegyHdrMap const* hdrMap );
  csSegyHdrMap();
  ~csSegyHdrMap();

  void readHdrMapExternalFile( std::string filename );

  /**
  * @return ID number of current default map
  */
  int mapID() const { return myDefaultMap; }
  /**
  * Add header to SEG-Y map
  * @param byteLoc   Byte location (0 for first byte)
  * @param byteSize  Header size in bytes
  * @param inType    Type of SEG-Y trace header as it shall be decoded from binary file
  * @param hdr       Header info defining output trace header
  */
  bool addHeader( int byteLoc, int byteSize, type_t inType, cseis_geolib::csHeaderInfo const& hdr );
  /**
  * Add header to SEG-Y map
  * @param hdr       SEG-Y header info defining SEG-Y trace header
  */
  bool addHeader( cseis_geolib::csSegyHeaderInfo const& hdr );
  /**
  * Add header to SEG-Y map
  * @param theName      Name of trace header
  * @param theByteLoc   Byte location (0 for first byte)
  * @param theByteSize  Header size in bytes
  * @param theInType    Type of SEG-Y trace header as it shall be decoded from binary file
  */
  bool addHeader( std::string const& theName, int theByteLoc, int theByteSize, type_t theInType );
  /**
  * Add header to SEG-Y map
  * @param theName      Name of trace header
  * @param theByteLoc   Byte location (0 for first byte)
  * @param theByteSize  Header size in bytes
  * @param theInType    Type of SEG-Y trace header as it shall be decoded from binary file
  * @param theDesc      Description of trace header
  */
  bool addHeader( std::string const& theName, int theByteLoc, int theByteSize, type_t theInType, std::string const& theDesc );
  /**
  * Add header to SEG-Y map
  * @param theName      Name of trace header
  * @param theByteLoc   Byte location (0 for first byte)
  * @param theByteSize  Header size in bytes
  * @param theInType    Type of SEG-Y trace header as it shall be decoded from binary file
  * @param theOutType   Type of output trace header
  * @param theDesc      Description of trace header
  */
  bool addHeader( std::string const& theName, int theByteLoc, int theByteSize, type_t theInType, type_t theOutType, std::string const& theDesc );
  /**
  * Add header to SEG-Y map
  * @return true if header was successfully removed
  */
  bool addHeader( int byteLoc, std::string const& typeNameIn, std::string const& typeNameOut, std::string const& name, std::string const& desc );

  bool removeHeader( int hdrIndex );
  bool removeHeader( std::string );
  void removeAll();
  int  headerIndex( std::string hdrName ) const;
  bool contains( std::string hdrName, int* hdrIndex ) const;
  /**
  * @param hdrIndex     Index of header in SEG-Y header map, starting with 0.
  * @return SEG-Y header definition
  */
  cseis_geolib::csSegyHeaderInfo const* header( int hdrIndex ) const;
  cseis_geolib::csSegyHeaderInfo const* header( std::string hdrName ) const;
  /**
  * @return Number of defined trace headers
  */
  int numHeaders() const;

  // Special methods for application of coordinate scalars
  int numCoordHeaders() const;
  int numElevHeaders() const;
  int numStatHeaders() const;

  int getCoordScalarHeaderIndex() const;
  int getElevScalarHeaderIndex() const;
  int getStatScalarHeaderIndex() const;

  int getCoordHeaderIndex( int index ) const;
  int getElevHeaderIndex( int index ) const;
  int getStatHeaderIndex( int index ) const;

  /**
  * Initialization of coordinate and elevation scalars
  * (1) Determine receiver & shot coordinates and elevations
  * (2) Determine whether SEG-Y coordinate and/or elevation scalars are specified in SEG-Y header map.
  * Call this function if trace headers concerning elevation, coordinate and/or scalars were changed,
  * for example by manually setting non-standard headers
  */
  void initScalars();
  
  void applyCoordinateScalar( cseis_geolib::csFlexHeader* hdrValues ) const;
//  void applyCoordinateScalar( cseis_geolib::csFlexHeader* hdrValues, int scalarCoord, int scalarElev ) const;
  void applyCoordinateScalarWriting( cseis_geolib::csFlexHeader* hdrValues ) const;

  void dump( FILE* stream ) const;

private:
  void init( bool initScalars );
  bool getHdrIndex( std::string const& name, int* hdrIndex ) const;
  bool getHdrIndex( std::string const& name, int* hdrIndex, int byteLocIn, int* hdrIndexByteLoc ) const;

  int myDefaultMap;
  bool myIsInitScalars;
  bool myReplaceExistingHeader;
  cseis_geolib::csVector<cseis_geolib::csSegyHeaderInfo const*>* myHdrList;

  cseis_geolib::csVector<int>* myElevHdrIDList;
  cseis_geolib::csVector<int>* myCoordHdrIDList;
  cseis_geolib::csVector<int>* myStatHdrIDList;
  int myHdrID_scalar_elev;
  int myHdrID_scalar_coord;
  int myHdrID_scalar_stat;
};

} // end namespace
#endif

