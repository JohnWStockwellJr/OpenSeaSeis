/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSegdHeader_DIGISTREAMER.h"
#include "geolib_endian.h"
#include <ostream>
#include <cstring>

using namespace std;
using namespace cseis_segd;

csExtendedHeader_DIGI::csExtendedHeader_DIGI() : csExtendedHeader() {
  myIsLittleEndian = cseis_geolib::isPlatformLittleEndian();
  numElements = 0;
  text = NULL;
}
csExtendedHeader_DIGI::~csExtendedHeader_DIGI() {
  if( text != NULL ) {
    delete [] text;
    text = NULL;
  }
}
void csExtendedHeader_DIGI::extractHeaders( byte const* buffer, int numBytes ) {
  numElements = numBytes;
  text = new char[numElements];
  memcpy(text,buffer,numElements);
}
//----------------------------------------------------------------------

csTraceHeaderExtension_DIGI::csTraceHeaderExtension_DIGI( int numBlocks ) : csTraceHeaderExtension( numBlocks ) {
  myIsLittleEndian = cseis_geolib::isPlatformLittleEndian();
}
csTraceHeaderExtension_DIGI::~csTraceHeaderExtension_DIGI() {
}
void csTraceHeaderExtension_DIGI::extractHeaders( byte const* buffer, commonTraceHeaderStruct* comTrcHdr ) {
  int blockSize = 32;

  comTrcHdr->rcvLineNumber  = UINT24(&buffer[0]);
  comTrcHdr->rcvPointNumber = UINT24(&buffer[3]);
  comTrcHdr->rcvPointIndex  = UINT8(&buffer[6]);

  if( myNumBlocks < 2 ) return;
  if( myIsLittleEndian ) {
    char buffer8[8];
    memcpy( buffer8, &buffer[0+blockSize], 8 );
    cseis_geolib::swapEndian( buffer8, 8, 8 );
    memcpy( &comTrcHdr->rcvEasting, buffer8, 8 );

    memcpy( buffer8, &buffer[8+blockSize], 8 );
    cseis_geolib::swapEndian( buffer8, 8, 8 );
    memcpy( &comTrcHdr->rcvNorthing, buffer8, 8 );

    memcpy( &buffer8, &buffer[16+blockSize], 4 );
    cseis_geolib::swapEndian4( buffer8, 4 );
    memcpy( &comTrcHdr->rcvElevation, buffer8, 4 );
  }
  else {
    memcpy( &comTrcHdr->rcvEasting, &buffer[0+blockSize], 8 );
    memcpy( &comTrcHdr->rcvNorthing, &buffer[8+blockSize], 8 );
    memcpy( &comTrcHdr->rcvElevation, &buffer[16+blockSize], 4 );
  }

  if( myNumBlocks < 6 ) return;
  comTrcHdr->serialNumber = UINT24(&buffer[1+blockSize*5]);
}
//--------------------------------------------------------
//
void csTraceHeaderExtension_DIGI::extractHeaders( byte const* buffer ) {
  memset( &block2, 0, sizeof(su_header) );

  block1.rcvLineNumber  = UINT24(&buffer[0]);
  block1.rcvPointNumber = UINT24(&buffer[3]);
  block1.rcvPointIndex  = UINT8(&buffer[6]);
  block1.numSamples     = UINT24(&buffer[7]);
  block1.extendedRcvLineNumber  = UINT32(&buffer[10]);
  block1.extendedRcvPointNumber = UINT32(&buffer[15]);
  block1.sensorType = UINT8(&buffer[20]);

  if( myNumBlocks < 2 ) return;
  if( myIsLittleEndian ) {
    cseis_geolib::swapEndian( (char*)&buffer[32], 2, 240 );
  }
  memcpy( &block2, &buffer[32], 240 );
  memcpy( block3.buffer, &buffer[272], 16 );
}

namespace cseis_segd {
void csExtendedHeader_DIGI::dump( std::ostream& os )
{
  os <<
    header("ExtendedHeader start") << '\n';
  if( text != NULL ) {
    os << text << endl;
  }
  os << header("ExtendedHeader end  ");
}

void csTraceHeaderExtension_DIGI::dump( std::ostream& os ) {
  os <<
    header("TraceHeaderExtension start") << endl;
  {
    os <<
      "rcvLineNumber          : " << block1.rcvLineNumber << endl <<
      "rcvPointNumber         : " << block1.rcvPointNumber << endl <<
      "rcvPointIndex          : " << block1.rcvPointIndex << endl <<
      "numSamples             : " << block1.numSamples << endl <<
      "extendedRcvLineNumber  : " << block1.extendedRcvLineNumber << endl <<
      "extendedRcvPointNumber : " << block1.extendedRcvPointNumber << endl <<
      "sensorType             : " << block1.sensorType << endl;
    if( myNumBlocks < 2 ) goto endofdump;
    os << " SU header: " << endl <<
      " Trace sequence number within line:        " << block2.tracl << endl << 
      " Trace sequence number within SEG Y file:  " << block2.tracr << endl <<
      " Original field record number:             " << block2.fldr << endl <<
      " Trace number within original field record:" << block2.tracf << endl <<
      " Energy source point number:               " << block2.ep    << endl <<
      " Ensemble number (i.e. CDP, CMP, CRP,...): " << block2.cdp   << endl <<
      " Trace number within the ensemble:         " << block2.cdpt  << endl <<
      " Trace identification code:                " << block2.trid  << endl <<
      " Number of vertically summed traces:       " << block2.nvs   << endl <<
      " Number of horizontally summed traces:     " << block2.nhs   << endl <<
      " Data use:                                 " << block2.duse  << endl <<
      " Distance from the center of the source point: " << block2.offset << endl <<
      " Receiver group elevation from sea level:  " << block2.gelev  << endl <<
      " Surface elevation at source:              " << block2.selev  << endl <<
      " Source depth below surface:               " << block2.sdepth << endl <<
      " Datum elevation at receiver group:        " << block2.gdel   << endl <<
      " Datum elevation at source:                " << block2.sdel   << endl <<
      " Water depth at source:                    " << block2.swdep  << endl <<
      " Water depth at receiver group:            " << block2.gwdep  << endl <<
      " Scalar to be applied to the previous 7 entries: " << block2.scalel << endl <<
      " Scalar to be applied to the next 4 entries:     " << block2.scalco << endl <<
      " Source coordinate - X:                    " << block2.sx    << endl <<
      " Source coordinate - Y:                    " << block2.sy    << endl <<
      " Group coordinate - X:                     " << block2.gx    << endl <<
      " Group coordinate - Y:                     " << block2.gy    << endl <<
      " Coordinate units:                         " << block2.counit << endl <<

      " mute time--start:                         " << block2.muts << endl <<
      " mute time--end:                           " << block2.mute << endl << 
      " number of samples in this trace :         " << block2.ns << endl <<
      " sample interval; in micro-seconds:        " << block2.dt << endl <<
      " gain type of field instruments code:      " << block2.gain << endl <<
      " instrument gain constant:                 " << block2.igc << endl <<
      " instrument early or initial gain:         " << block2.igi << endl <<
      " alias filter frequency if used:           " << block2.afilf << endl <<
      " alias filter slope:                       " << block2.afils << endl <<
      " notch filter frequency if used:           " << block2.nofilf << endl <<
      " notch filter slope:                       " << block2.nofils << endl <<
      " low cut frequency if used:                " << block2.lcf << endl <<
      " high cut frequncy if used:                " << block2.hcf << endl <<
      " low cut slope:                            " << block2.lcs << endl <<
      " high cut slope:                           " << block2.hcs << endl <<
      " year data recorded:                       " << block2.year  << endl <<
      " day of year:                              " << block2.day   << endl <<
      " hour of day (24 hour clock):              " << block2.hour  << endl <<
      " minute of hour:                           " << block2.minute << endl <<
      " second of minute:                         " << block2.sec    << endl <<
      " time basis code:                          " << block2.timbas << endl <<
      " trace weighting factor, defined as 1/2^N: " << block2.trwf << endl <<
      " geophone group number of roll switch:     " << block2.grnors << endl <<
      " geophone group number of trace one:       " << block2.grnofr << endl <<
      " geophone group number of last trace:      " << block2.grnlof << endl <<
      " gap size (total number of groups dropped):" << block2.gaps  << endl <<
      " overtravel taper code:                    " << block2.otrav << endl <<
      " sample spacing for non-seismic data:      " << block2.d1 << endl << 
      " first sample location for non-seismic data: " << block2.f1 << endl <<
      " sample spacing between traces:            " << block2.d2 << endl <<
      " first trace location:                     " << block2.f2 << endl <<
      " negative of power used for range compression: " << block2.ungpow << endl <<
      " reciprocal of scaling factor to normalize:    " << block2.unscale << endl <<
      " number of traces:                         " << block2.ntr << endl <<
      " mark selected traces:                     " << block2.mark << endl <<
      endl;
  }
 endofdump:
  os << header("TraceHeaderExtension end  ") << endl;

}

}

