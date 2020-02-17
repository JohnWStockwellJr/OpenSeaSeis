/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TRACE_HEADER_DEF_H
#define CS_TRACE_HEADER_DEF_H

#include <string>
#include "geolib_defines.h"

namespace cseis_geolib {
  class csHeaderInfo;
  template<typename T> class csVector;
}

namespace cseis_system {

static int const HEADER_NOT_FOUND = -1;
static int const HEADER_EXISTS    = -2;
static int const HEADER_IS_VECTOR = -4;

//static int const HEADER_INDEX_ABSOLUTE_TIME_DAYS = 0;
//static int const HEADER_INDEX_ABSOLUTE_TIME_SEC  = 1;
//static int const HEADER_INDEX_ABSOLUTE_TIME_US   = 2;

  class csTraceHeaderInfo;
  class csTraceHeaderData;
  class csMemoryPoolManager;
  class csTraceGather;

/**
* Definition of trace header
*
* This class defines all headers in the trace header (as seen by one module)
* Each header is defined by a header info object:
* - Name
* - Type (int, float...)
* - Description
* - Number of elements (for string and array types)
*
* Each module has its own trace header definition object.
* Header info objects are only stored as const pointers. They are allocated and memory-managed in csTraceHeaderInfoPool
*
* @author Bjorn Olofsson
* @date   2007
*/
class csTraceHeaderDef {
  static int const HDR_ID_TIME_SAMP1_S  = 0;
  static int const HDR_ID_TIME_SAMP1_US = 1;

 public:
  /**
   * Use this constructor to create header definition for CSEIS flows
   * @param numInputPorts: Number of input ports
   * @param hdefPrev:   List of pointers to trace header definition objects of all input ports
   * @param memManager: Pointer to memory manager object
   */
  csTraceHeaderDef( int numInputPorts, csTraceHeaderDef const** hdefPrev, csMemoryPoolManager* memManager );
  /**
   * Use this constructor to create plain & empty header definition. This can be populated later using addHeader() functions
   */
  csTraceHeaderDef( csMemoryPoolManager* memManager );
  /**
   * Use this constructor to create header definitions for any other purpose
   * @param refHdrDefPtr:  Reference header definition. Required input parameter for internal settings.
   */
  csTraceHeaderDef( csTraceHeaderDef* refHdrDefPtr );
  ~csTraceHeaderDef();
  /**
  * Add new trace header
  * @param type  (i) Header type, e.g. TYPE_INT, TYPE_FLOAT, TYPE_DOUBLE, TYPE_STRING...
  * @param name  (i) Header name
  * @param nElements (i) Number of header elementsfor array types (for strign headers, supply maximum length of string here, must be larger than 1)
  * @return Access index of new trace header
  */
  int addHeader( cseis_geolib::type_t type, std::string const& name, int nElements = 1 );
  /**
  * Add new trace header, with description
  * @param type  (i) Header type, e.g. TYPE_INT, TYPE_FLOAT, TYPE_DOUBLE, TYPE_STRING...
  * @param name  (i) Header name
  * @param name  (i) Description
  * @param nElements (i) Number of header elementsfor array types (for strign headers, supply maximum length of string here, must be larger than 1)
  * @return Access index of new trace header
  */
  int addHeader( cseis_geolib::type_t type, std::string const& name, std::string const& description, int nElements = 1 );
  /**
  * Add new trace header
  * @param info  (i) Header info
  * @param nElements (i) Number of header elementsfor array types (for strign headers, supply maximum length of string here, must be larger than 1)
  * @return Access index of new trace header
  */
  int addHeader( cseis_geolib::csHeaderInfo const* info, int nElements = 1 );
  /**
  * Add standard trace header
  * @param name  (i) Header name
  * @return Access index of new trace header
  */
  int addStandardHeader( std::string const& name );
  /**
   * Delete trace header.
   * NOTE: Trace headers are not deleted for the current module. Instead they will be deleted for the following modules
   * @return old header index
   */
  int deleteHeader( std::string const& name );
  /**
   * Delete all trace headers
   * ..excluding required system headers
   */
  void deleteAllHeaders();
  /**
   * @return true if the two header definitions are equal, i.e. they contain the same trace header names with same sorting
  */
  bool equals( csTraceHeaderDef const* hdef ) const;
  /// Return true if header exists
  bool headerExists( std::string const& name ) const;
  /// Return sequential index of given header
  int headerIndex( std::string const& name ) const;
  /// Return type of given header (at given sequential index)
  cseis_geolib::type_t headerType( int index ) const;
  /// Return type of given header
  cseis_geolib::type_t headerType( std::string const& name ) const;
  csTraceHeaderInfo const* headerInfo( int index ) const;
  /// Return name of given header (at given sequential index)
  std::string headerName( int index ) const;
  /// Return description of given header (at given sequential index)
  std::string headerDesc( int index ) const;
  /// Return number of elements in this trace header (1 for normal headers, n>1 for array or string headers)
  int numElements( int index ) const;

  /// Return total number of headers
  int numHeaders() const;
  /// @return number of byte offset of specified header inside byte buffer
  int getByteLocation( int index ) const;
  /// Return total number of bytes required to store trace header values
  inline int getTotalNumBytes() const { return myTotalNumBytes; }
  /**
  * @return list of indices of trace headers that shall be deleted from the seismic traces that are input to a module
  */
  cseis_geolib::csVector<int> const* getIndexOfHeadersToDel() const;

  int mpi_getByteSize() const;
  int mpi_compress( char* data ) const;
  void mpi_decompress( char const* data );

  void dump( FILE* fout ) const;
  /**
  * Set byte location... Call before using this object in exec phase
  */
  void resetByteLocation();
  /**
  * Retrieve pointer handle to list of byte locations for each trace header value in trace header value block.
  */
  inline int const* getHandleByteLocation() const {
    return myByteLocation;
  }
  bool isSystemTraceHeader( std::string const& name ) const;
  friend class csTraceGather;
private:
  static bool isSystemTraceHeader( int index );
  int addHeader_internal( cseis_geolib::type_t type, std::string const& name, std::string const& description, int nElements );
  /**
   * Delete trace header.
   * NOTE: Trace headers are not deleted for the current module. Instead they will be deleted for the following modules
   */
  void deleteHeader( int index );
  csMemoryPoolManager* getMemoryManager();

  csTraceHeaderDef();
  csTraceHeaderDef( csTraceHeaderDef const& obj );
  csTraceHeaderDef& operator=( const csTraceHeaderDef& obj );
  
  /// Manages trace header info objects
  csMemoryPoolManager* myMemoryManager;
  /// List containing pointers to info objects defining all trace headers in this definition object
  cseis_geolib::csVector<csTraceHeaderInfo const*>* myTraceHeaderInfoList;

  /// List of sequential header indexes of all trace headers that need to be added (one list per input port)
  cseis_geolib::csVector<csTraceHeaderInfo const*>* myTraceHeadersToAdd;
  cseis_geolib::csVector<int>* myIndexOfHeadersToDel;

  /// Number of input ports. Each input port may have its own set of trace headers which need to be synchronized
  int myNumInputPorts;
  /// Total number of bytes required to store trace header values
  int myTotalNumBytes;

  /// TEST: Copy of myByteLocation array in csTraceHeaderData. Try to do the same thing in this class, replace the other if possible
  /// Maps the sequential header index 0,1,2... to the byte location index in the char* 'value block'
  int* myByteLocation;

  /// Return index to given header. The index can be used in all setter and getter methods in the trace data object (defined separately)
  bool getIndex( std::string const& name, int& index ) const;
  /// Initialize object with header definitions from all input ports
  /// @param hdefPrev:   List of pointers to trace header definition objects of all input ports
  void initInputPorts( csTraceHeaderDef const** hdefPrev );
  /// General initialisation
  void init( int numInputPorts, csMemoryPoolManager* memManager );
  static std::string HEADER_NAME_UNKNOWN;

};

} // namespace
#endif

