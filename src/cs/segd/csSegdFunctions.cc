/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSegdDefines.h"
#include "csSegdFunctions.h"
#include <iomanip>
#include <cmath>

using namespace std;

namespace cseis_segd {
  int bcd( byte const* bytes, int startNibble, int numNibbles )
  {
    // One nibble consists of 4 bits = half a byte
    int val = 0;
    int startByte = (int)(startNibble / 2);
    int endNibblePlus = startNibble+numNibbles;

    for( int iNibble = startNibble, iByte = startByte; iNibble < endNibblePlus; iNibble++ ) {
      val *= 10;
      if( iNibble & 0x1 ) {
        val += bytes[iByte++] & 0xf;
        //      fprintf(stderr,"A: #%d   #%d :  %x  %d\n", iNibble, iByte, bytes[iByte] & 0xf, bytes[iByte] & 0xf);
      }
      else {
        val += (bytes[iByte] >> 4) & 0xf;
        //      fprintf(stderr,"A: #%d   #%d :  %x  %d\n", iNibble, iByte, (bytes[iByte]>>4) & 0xf, (bytes[iByte] >> 4 ) & 0xf);
      }
    }
    return val;
  }

  //---------------------------------------------------------------

  int string2int( byte const* str, int length ) {
    return (int)string2long( str, length );
  }

  long string2long( byte const* str, int length ) {
    int digit;
    long result = 0;
    bool isNegative = false;

    //  int length = str.size();
    for( int i = 0; i < length; i++ ) {
      if( str[i] >= '0' && str[i] <= '9' ) {
        digit = str[i] - '0';
        result *= 10;
        result += digit;
      }
      else if( str[i] == '-' ) {
        if( result != 0 ) return 0; // Wrong format!
        isNegative = true;
      }
      else if( str[i] == '.' ) {
        break;
      }
      else if( str[i] != ' ' ) {
        return 0;
      }
    }
    if( isNegative ) {
      result = -result;
    }
    return result;
  }

  //-------------------------------------------------------------
  //
  double string2double(byte const* buffer, int length) {
    double result = 0.0;
    int digit;
    bool isNegative = false;

    int counter = 0;
    // Remove leading blanks
    while( counter < length && buffer[counter] == ' ' ) {
      counter++;
    }
    if( counter == length ) return result;
  
    // Check for sign
    if( buffer[counter] == '-' ) {
      isNegative = true;
      counter++;
    }
  
    // Leading digits
    while( counter < length && (buffer[counter] >= '0' && buffer[counter] <= '9') ) {
      digit = buffer[counter] - '0';
      result *= 10.0;
      result += digit;
      counter++;
    }
  
    // Optional dot for floating point number
    if( counter < length && buffer[counter] == '.' ) {
      counter++;
      // Trailing digits
      double decimalMultiplier = 0.1;
      while( counter < length && (buffer[counter] >= '0' && buffer[counter] <= '9') ) {
        digit = buffer[counter] - '0';
        result += decimalMultiplier*(double)digit;
        decimalMultiplier *= 0.1;
        counter++;
      }
    }
  
    // Optional 'e' notation
    if( counter < length && (buffer[counter] == 'e' || buffer[counter] == 'E') ) {
      counter++;
      bool isExponentNegative = false;
      double exponent = 0.0;
    
      if( counter >= length || ( buffer[counter] != '+' && buffer[counter] != '-' ) ) {
        return 0.0;  // Format error
      }
      if( buffer[counter] != '-' ) {
        isExponentNegative = true;
      }
      counter++;
      int checkIndex = counter;
      while( counter < length && (buffer[counter] >= '0' && buffer[counter] <= '9') ) {
        digit = buffer[counter] - '0';
        exponent *= 10.0;
        exponent += digit;
        counter++;
      }
      if( counter == checkIndex ) {
        //throw( ExpressionException( "Wrong number format (no number in 'e' notation)",
        return 0.0;
      }
      if( isExponentNegative ) exponent = -exponent;
      result = pow( result, exponent );
    }
    if( counter != length ) {  // Format error
      return 0.0;
    }
    if( isNegative ) result = -result;
    return result;
  }

  float string2float(byte const* buffer, int length) {
    return (float)string2double( buffer, length );
  }

  csInt64_t convert2uint64( byte const* buffer ) {
    return( ((csInt64_t)*(buffer+7)) +
            ((csInt64_t)*(buffer+6)<<8UL) +
            ((csInt64_t)*(buffer+5)<<16UL) +
            ((csInt64_t)*(buffer+4)<<24UL) +
            ((csInt64_t)*(buffer+3)<<32UL) +
            ((csInt64_t)*(buffer+2)<<40UL) +
            ((csInt64_t)*(buffer+1)<<48UL) +
            ((csInt64_t)*(buffer+0)<<56UL ) );
  }

  //--------------------------------------------------------------
  int sampleBitSize( int segdFormatCode ) {
    switch( segdFormatCode ) {
    case 8015: // 20 bit binary demultiplexed
      return 20;
    case 8022: //       8        bit quaternary demultiplexed
      return 8;
    case 8024: //               16 bit quaternary demultiplexed
      return 16;
    case 8036: //               24 bit 2's compliment integer demultiplexed
      return 24;
    case 8038: //               32 bit 2's compliment integer demultiplexed
      return 32;
    case 8042: //               8 bit hexadecimal demultiplexed
      return 8;
    case 8044: //               16 bit hexadecimal &multiplexed
      return 16;
    case 8048: //               32 bit hexadecimal demultiplexed
      return 32;
    case 8058: //               32 bit IEEE demultiplexed
      return 32;
    default:
      return 0; // not supported
    }
  }
  bool isFormatCodeSupported( int segdFormatCode ) {
    switch( segdFormatCode ) {
    case 8022: //       8        bit quaternary demultiplexed
    case 8024: //               16 bit quaternary demultiplexed
    case 8038: //               32 bit 2's compliment integer demultiplexed
    case 8042: //               8 bit hexadecimal demultiplexed
    case 8044: //               16 bit hexadecimal &multiplexed
    case 8048: //               32 bit hexadecimal demultiplexed
      return false;
    case 8015: // 20 bit binary demultiplexed
    case 8036: // 24 bit 2's compliment integer demultiplexed
    case 8058: // 32 bit IEEE demultiplexed
      return true;
    default:
      return false;
    }
  }
  bool isManufacturerSupported( int manufactCode ) {
    switch( manufactCode ) {
    case 9:  // geospace
    case 13:  // Sercel
    case 18:  // I/O
      return true;
    case 20:  // Fairfield
      return true;
    case 22:  // Geco-Prakla
      return false;
    default:
      return false;
    }
  }
  int manufacturerRecordingSystem( int manufactCode ) {
    switch( manufactCode ) {
    case 9:  // geospace
      return RECORDING_SYSTEM_GEORES;
    case 13:  // Sercel
      return RECORDING_SYSTEM_SEAL;
    case 18:  // I/O
      return RECORDING_SYSTEM_DIGISTREAMER;
    case 20:  // Fairfield
      return RECORDING_SYSTEM_FAIRFIELD;
    case 22:  // Geco Prakla
      return RECORDING_SYSTEM_GECO;
    case 42: // Gunlink SEGD?
      return UNKNOWN;
    default:
      return UNKNOWN;
    }
  }
  void getManufacturerName( int manufactCode, string* text ) {
    switch( manufactCode ) {
    case 2:
      *text = "Applied Magnetics Corporation";
      break;
    case 3:
      *text = "Litton Resources Systems, Inc.";
      break;
    case 4:
      *text = "SIE, Inc.";
      break;
    case 7:
      *text = "Electro-Technical Labs";
      break;
    case 9:
      *text = "Geo Space Corporation ";
      break;
    case 13:
      *text = "Sercel (Societe d'Etudes, Recherches Et Constructions Electroniques)";
      break;
    case 15:
      *text = "Texas Instruments, Inc.";
      break;
    case 17:
      *text = "GUS Manufacturing, Inc.";
      break;
    case 18:
      *text = "Input/Output, Inc.";
      break;
    case 19:
      *text = "Geco-Prakla";
      break;
    case 20:
      *text = "Fairfield";
      break;
    case 22:
      *text = "Geco-Prakla";
      break;
    case 32:
      *text = "Halliburton";
      break;
    case 33:
      *text = "CompuSeis";
      break;
    case 34:
      *text = "Syntron Inc";
      break;
    case 35:
      *text = "Syntron Europe Ltd";
      break;
    case 39:
      *text = "Grant Geophysical";
      break;
    default:
      *text = "UNKNOWN";
    }
  }
  bool isRevisionSupported( int rev1, int rev2 ) {
    if( ((rev1 == 0 || rev1 == 1 || rev1 == 2) && rev2 == 0 ) ||
        (rev1 == 2 && rev2 == 1) || (rev1 == 0 && rev2 == 1) ) {
      return true;
    }
    else if( (rev1 == 0 || (rev1==1 && (rev2 == 5 || rev2 == 6)) )) { // Rev 1.5 & 1.6 (Fairfield)
      return true;
    }
    else {
      return false;
    }
  }

  void dumpRawHex( std::ostream& os, byte const* buffer, int numBytes )
  {
    for( int i = 0; i < numBytes; i++ ) {
      //std::cout << int( (h.myBuffer[i] & 0xf0) >> 4 ) << int(h.myBuffer[i] & 0x0f) << ' ';
      os << std::hex << std::setw(2) << std::setfill('0') << int(buffer[i]) << ' ' ;
      if (i != 0 && (i+1) % (numBytes/2) == 0)
        os << std::endl;
    }
    os << std::dec;
  }

  void dumpRawHex( FILE* file, byte const* buffer, int numBytes )
  {
    for( int i = 0; i < numBytes; i++ ) {
      fprintf(file,"%02x ", (int)buffer[i] );
      if( i != 0 && (int)((i+1)/16)*16 == (i+1) ) { 
        fprintf(file,"\n");
      }
    }
  }

  void dumpRawASCII( std::ostream& os, byte const* buffer, int numBytes )
  {
    int num = int((numBytes-1)/120) + 1;
    int counterTotal = 0;
    for( int i = 0; i < num; i++ ) {
      int counter = 0;
      while( counterTotal < numBytes && counter < 120 ) {
        os << (char)buffer[counterTotal];
        counter++;
        counterTotal++;
      }
      os << endl;
    }
  }

  void dumpCommonHeaders( commonTraceHeaderStruct& comTrcHdr ) {
    fprintf(stdout,"chanNum: %d, chanTypeID: %d, traceEdit: %d\n", comTrcHdr.chanNum, comTrcHdr.chanTypeID, comTrcHdr.traceEdit);
    fprintf(stdout,"rcvLineNum: %d, rcvPointNum: %d, rcvPointIndex: %d\n", comTrcHdr.rcvLineNumber, comTrcHdr.rcvPointNumber, comTrcHdr.rcvPointIndex );
    fprintf(stdout,"rcvEast: %f, rcvNorth: %f, rcvElev: %f, serialNum: %d\n", comTrcHdr.rcvEasting, comTrcHdr.rcvNorthing, comTrcHdr.rcvElevation, comTrcHdr.serialNumber );
  }

} // end namespace

