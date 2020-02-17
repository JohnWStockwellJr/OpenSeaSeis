/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_IO_DEFINES_H
#define CS_IO_DEFINES_H

namespace cseis_io {

  typedef unsigned char byte;
  static const char ID_TEXT_CSEIS[] = "CSEIS";
  static const char ID_TEXT_OSEIS[] = "OSEIS";

 class csIODefines {
 public:
   static void createVersionString( int versionNumber, char* versionText4Chars ) {
     versionNumber = versionNumber % 100;
     int versionMajor = (int)(versionNumber/10);
     int versionMinor = (int)(versionNumber - versionMajor*10);
     sprintf( versionText4Chars, "%1d.%1d", versionMajor, versionMinor );
     versionText4Chars[3] = '\0';
   }
 };

} // end namespace
#endif

