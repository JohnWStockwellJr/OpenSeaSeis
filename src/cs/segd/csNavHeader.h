/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_NAV_HEADER_H
#define CS_NAV_HEADER_H

#include "csSegdDefines.h"
#include <ostream>

namespace cseis_segd {

struct hdrDesc_HYDRONAV_Struct {
  int version;         // 1
  int revision;        // 2
  int length;          // 3-4
  int numTagBlocks;    // 5
  int quality;         // 6
};

struct tagDesc_HYDRONAV_Struct {
  char desc[8+1];   //  1-8  Tag description */
  int  version;     //  9-9  Tag data version */
  int  revision;    // 10-10 Tag data revision */
  int  length;      // 11-12 Length of tag block */
  int  dateFormat;  // 13-13 Date format. if tag == HYDRONAV: 0x00 means MM/DD/YYYY, 0x01 means DD/MM/YYYY
};

typedef struct
{
  hdrDesc_HYDRONAV_Struct hdrDesc;    //   hdrDesc
  tagDesc_HYDRONAV_Struct tagDesc;    //   tagDesc
  char header[8+1];        //  17-24  Hydronav header "@@@UHEAD"
  char time[15+1];         //  25-39  Time in HH:MM:SS:SSSSSS (A15)
  int  GPSclock;           //  41-41  GPS clock (I1)
  char vesselId[8+1];      //  43-50 Vessel id (A8)
  int  stationNum;         //  52-58  Station number (I7) 
  int  FSID;               //  60-69  FSID (I10) 
  double northing;         //  71-82  Source northing (F11.2N) 
  double easting;          //  84-94  Source easting (F11.2E) 
  float vesselSpeed;       //  97-102 Vessel speed (F6.2) 
  float bearing;           // 106-111 bearing (F6.2) 
  float waterDepth;        // 115-119 Water depth (F5.1) 
  char date[10+1];         // 122-131 Date (A10) 
  int distOffline;         // 133-138 Distance offline (I5M) 
  int distInline;          // 140-148 Distance inline (I8M) 
  int lineStatus;          // 150-150 Line status (I2)
  //  0 Off line
  //  1 Approach
  //  2 On line
  //  3 Run out 
  float GPSHeight;         // 152-157 GPS height (F5.1M) 
  char lineName[16+1];          //  159-174  Line name (A16) 
  char eod;                    // 175-175 End of HydroNav data '*' 
  char eot;                    // 176-176 End of tag mark 0xff 
  tagDesc_HYDRONAV_Struct tagDescGCS90;
}
NavHeader_HYDRONAV_ver6_Struct;  // Revision 6

typedef struct
{
  hdrDesc_HYDRONAV_Struct hdrDesc;    //   hdrDesc
  tagDesc_HYDRONAV_Struct tagDesc;    //   tagDesc
  char header[8+1];            //  17-24  Hydronav header "@@@UHEAD"
  char time[15+1];             //  25-39  Time in HH:MM:SS:SSSSSS (A15)
  int  GPSclock;               //  41-41  GPS clock (I1)
  char lineName[8+1];          //  43-50  Line name (A8) 
  int  stationNum;        //  52-58  Station number (I7) 
  int  FSID;             //  60-69  FSID (I10) 
  double northing;         //  71-81  Source northing (F11.2N) 
  double easting;          //  84-94  Source easting (F11.2E) 
  float vesselSpeed;       //  97-102 Vessel speed (F6.2) 
  float bearing;           // 106-111 bearing (F6.2) 
  float waterDepth;        // 115-119 Water depth (F5.1) 
  char date[10+1];             // 122-131 Date (A10) 
  int distOffline;       // 133-138 Distance offline (I5M) 
  int distInline;        // 140-148 Distance inline (I8M) 
  int lineStatus;        // 150-150 Line status (I2)
  //    0 Off line
  //    1 Approach
  //    2 On line
  //    3 Run out 
  float GPSHeight;         // 152-157 GPS height (F5.1M) 
  char vesselId[15+1];         // 159-173 Vessel id (A15) 
  char eod;                    // 175-175 End of HydroNav data '*' 
  char eot;                    // 176-176 End of tag mark 0xff 
}
NavHeader_HYDRONAV_ver1_Struct;

typedef struct {
  char   headerID[2+1];
  int    length;
  char   programRevision[4+1];
  int    lineStatus;
  char   shotTime[21+1];  // 2+2+2+1+6+4+2+2
  char   timeRef[3+1];
  int    shotNumber;
  char   lineName[16+1];
  double masterLatitude;
  double masterLongitude;
  float  waterDepth;
  double sourceLatitude;
  double sourceLongitude;
  float  masterGyro;
  float  masterCMG;
  float  masterSpeed;
}
NavHeader_PSI_Struct;

typedef struct {
  NavHeader_PSI_Struct   psi;
  double   masterEasting;
  double   masterNorthing;
  double   sourceDeltaEasting;
  double   sourceDeltaNorthing;
  float    lineBearing;
  int      julianDay;  
}
NavHeader_LABO_Struct;

/**
 * Segd Navigation header
 *
 * @author Bjorn Olofsson
 * @date 2006
 */
class csNavHeader {
public:
  csNavHeader( int navHeaderID );
  ~csNavHeader();
  void dump( std::ostream& os );
  void extractHeaders( byte const* buffer, int numBytes );
  void setCommonRecordHeaders( commonRecordHeaderStruct& comRecHdr );
  int size() const { return mySize; }
  int navHeaderID() const { return myNavHeaderID; }

private:
  csNavHeader( csNavHeader const& obj );
  void* myStructPtr;
  int myNavHeaderID;
  int mySize;
};

} // namespace
#endif

