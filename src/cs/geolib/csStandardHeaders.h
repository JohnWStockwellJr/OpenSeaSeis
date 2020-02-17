/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_STANDARD_HEADERS_H
#define CS_STANDARD_HEADERS_H

#include "geolib_defines.h"
#include "csHeaderInfo.h"

namespace cseis_geolib {

static csHeaderInfo const HDR_AN_AZIM       ( TYPE_FLOAT, "an_azim", "Azimuth angle from North [deg]" );
static csHeaderInfo const HDR_AN_INCI       ( TYPE_FLOAT, "an_inci", "Incidence angle from vertical [deg]" );
static csHeaderInfo const HDR_AN_ROLL       ( TYPE_FLOAT, "an_roll", "Roll angle [deg]" );
static csHeaderInfo const HDR_AN_TILT       ( TYPE_FLOAT, "an_tilt", "Tilt angle [deg]" );
static csHeaderInfo const HDR_AN_TILTX      ( TYPE_FLOAT, "an_tiltx", "Tilt angle (X direction) [deg]" );
static csHeaderInfo const HDR_AN_TILTY      ( TYPE_FLOAT, "an_tilty", "Tilt angle (Y direction) [deg]" );
static csHeaderInfo const HDR_AOFFSET       ( TYPE_FLOAT, "aoffset", "Absolute source-receiver offset [m]" );

static csHeaderInfo const HDR_BIN_X         ( TYPE_DOUBLE, "bin_x", "Bin X coordinate [m]" );
static csHeaderInfo const HDR_BIN_Y         ( TYPE_DOUBLE, "bin_y", "Bin Y coordinate [m]" );

static csHeaderInfo const HDR_CBL_AZIM      ( TYPE_FLOAT, "cbl_azim", "Cable azimuth form North [deg]");
static csHeaderInfo const HDR_CCP           ( TYPE_INT, "ccp", "Common conversion number" );
static csHeaderInfo const HDR_CHAN          ( TYPE_INT, "chan", "Channel number" );
static csHeaderInfo const HDR_CHAN_ENS      ( TYPE_INT, "chan_ens", "Trace/channel number within ensemble" );
static csHeaderInfo const HDR_CMP           ( TYPE_INT, "cmp", "CMP number" );
static csHeaderInfo const HDR_CMP_X         ( TYPE_DOUBLE, "cmp_x", "CMP X coordinate [m]" );
static csHeaderInfo const HDR_CMP_Y         ( TYPE_DOUBLE, "cmp_y", "CMP Y coordinate [m]" );
static csHeaderInfo const HDR_CMP_Z         ( TYPE_FLOAT, "cmp_z", "CMP Z coordinate [m]" );
static csHeaderInfo const HDR_COL           ( TYPE_INT, "col", "Column/crossline number" );

static csHeaderInfo const HDR_DATA_TYPE     ( TYPE_INT, "data_type", "Data type" );
static csHeaderInfo const HDR_DELAY_TIME    ( TYPE_FLOAT, "delay_time", "Delay time [ms]" );
static csHeaderInfo const HDR_DC            ( TYPE_FLOAT, "dc", "DC bias" );

static csHeaderInfo const HDR_FFID          ( TYPE_INT, "ffid", "Field file ID number" );
static csHeaderInfo const HDR_FILENO        ( TYPE_INT, "fileno", "Sequential input file number" );
static csHeaderInfo const HDR_FILT_HIGH_FREQ( TYPE_INT, "filt_high_freq", "High-cut filter [Hz]" );
static csHeaderInfo const HDR_FILT_HIGH_DB  ( TYPE_INT, "filt_high_db", "High-cut filter [dB/oct]" );
static csHeaderInfo const HDR_FILT_LOW_FREQ ( TYPE_INT, "filt_low_freq", "Low-cut filter [Hz]" );
static csHeaderInfo const HDR_FILT_LOW_DB   ( TYPE_INT, "filt_low_db", "Low-cut filter [dB/oct]" );
static csHeaderInfo const HDR_FOLD          ( TYPE_INT, "fold", "Trace fold" );
static csHeaderInfo const HDR_FOLD_VERT     ( TYPE_INT, "fold_vert", "Trace fold (vertically stacked)" );

static csHeaderInfo const HDR_GAIN_TYPE     ( TYPE_INT, "gain_type", "Gain type (1=fixed, 2=binary, 3=float, 4-N=optional)");
static csHeaderInfo const HDR_GAIN          ( TYPE_INT, "gain", "Gain constant [dB]");
static csHeaderInfo const HDR_GUN_SEQ       ( TYPE_INT, "gun_seq", "Gun sequence number" );

static csHeaderInfo const HDR_HEADING       ( TYPE_FLOAT, "heading", "Heading [deg]" );

static csHeaderInfo const HDR_INCL_I        ( TYPE_FLOAT, "incl_i", "Inclinometer inline component" );
static csHeaderInfo const HDR_INCL_C        ( TYPE_FLOAT, "incl_c", "Inclinometer crossline component" );
static csHeaderInfo const HDR_INCL_V        ( TYPE_FLOAT, "incl_v", "Inclinometer vertical component" );

static csHeaderInfo const HDR_MUTE_START    ( TYPE_FLOAT, "mute_start", "Start time of mute [ms]" );
static csHeaderInfo const HDR_MUTE_END      ( TYPE_FLOAT, "mute_end", "End time of mute [ms]" );

static csHeaderInfo const HDR_NODE          ( TYPE_INT, "node", "Node ID" );
static csHeaderInfo const HDR_NSAMP         ( TYPE_INT, "nsamp", "Number of live samples" );

static csHeaderInfo const HDR_OFFSET       ( TYPE_FLOAT, "offset", "Source-receiver offset (signed) [m]" );
static csHeaderInfo const HDR_ORIENT_I     ( TYPE_FLOAT, "orient_i", "Sensor orientation 'I' value" );
static csHeaderInfo const HDR_ORIENT_C     ( TYPE_FLOAT, "orient_c", "Sensor orientation 'C' value" );
static csHeaderInfo const HDR_ORIENT_V     ( TYPE_FLOAT, "orient_v", "Sensor orientation 'V' value" );

static csHeaderInfo const HDR_RCV          ( TYPE_INT,    "rcv", "Receiver station number" );
static csHeaderInfo const HDR_REC_INDEX    ( TYPE_INT,    "rec_index", "Receiver station index" );   // BUGFIX 21 May 2009
static csHeaderInfo const HDR_REC_LINE     ( TYPE_INT,    "rec_line", "Receiver line number" );
static csHeaderInfo const HDR_REC_X        ( TYPE_DOUBLE, "rec_x", "Receiver X coordinate [m]" );
static csHeaderInfo const HDR_REC_Y        ( TYPE_DOUBLE, "rec_y", "Receiver Y coordinate [m]" );
static csHeaderInfo const HDR_REC_Z        ( TYPE_FLOAT,  "rec_z", "Receiver Z coordinate [m]" );
static csHeaderInfo const HDR_REC_DATUM    ( TYPE_FLOAT,  "rec_datum", "Datum elevation/tide at receiver [m]" );
static csHeaderInfo const HDR_REC_ELEV     ( TYPE_FLOAT,  "rec_elev", "Receiver elevation [m]" );
static csHeaderInfo const HDR_REC_WDEP     ( TYPE_FLOAT,  "rec_wdep", "Water depth at receiver [m]" );
static csHeaderInfo const HDR_ROW           ( TYPE_INT, "row", "Row/inline number" );

static csHeaderInfo const HDR_SAMPINT_US   ( TYPE_INT, "sampint_us", "Sample interval [us]" );
static csHeaderInfo const HDR_SCALAR_ELEV  ( TYPE_FLOAT, "scalar_elev", "Elevation scalar [SEGY]" );
static csHeaderInfo const HDR_SCALAR_COORD ( TYPE_FLOAT, "scalar_coord", "Coordinate scalar [SEGY]" );
static csHeaderInfo const HDR_SCALAR_STAT  ( TYPE_FLOAT, "scalar_stat", "Static header scalar [SEGY]" );
static csHeaderInfo const HDR_SENSOR       ( TYPE_INT, "sensor", "Sensor number" );
static csHeaderInfo const HDR_SEQ          ( TYPE_INT, "seq", "Sequence number" );
static csHeaderInfo const HDR_SERIAL       ( TYPE_INT, "serial", "Sensor serial number" );
static csHeaderInfo const HDR_SAIL_LINE    ( TYPE_INT, "sail_line", "Sail line number" );
static csHeaderInfo const HDR_SOURCE       ( TYPE_INT, "source", "Source station/Shot point number" );
static csHeaderInfo const HDR_SOU_INDEX    ( TYPE_INT, "sou_index", "Source station/Shot point index" );
static csHeaderInfo const HDR_SOU_LINE     ( TYPE_INT, "sou_line", "Source line number" );
static csHeaderInfo const HDR_SOU_X        ( TYPE_DOUBLE, "sou_x", "Source X coordinate [m]" );
static csHeaderInfo const HDR_SOU_Y        ( TYPE_DOUBLE, "sou_y", "Source Y coordinate [m]" );
static csHeaderInfo const HDR_SOU_Z        ( TYPE_FLOAT, "sou_z", "Source Z coordinate [m]" );
static csHeaderInfo const HDR_SOU_DATUM    ( TYPE_FLOAT, "sou_datum", "Datum elevation/tide at source [m]" );
static csHeaderInfo const HDR_SOU_ELEV     ( TYPE_FLOAT, "sou_elev", "Source elevation [m]" );
static csHeaderInfo const HDR_SOU_WDEP     ( TYPE_FLOAT, "sou_wdep", "Water depth at source [m]" );
static csHeaderInfo const HDR_SR_AZIM      ( TYPE_FLOAT, "sr_azim", "Source-receiver azimuth from North [deg]" );

static csHeaderInfo const HDR_STAT_REC     ( TYPE_FLOAT, "stat_rec", "Receiver statics [ms]" );
static csHeaderInfo const HDR_STAT_SOU     ( TYPE_FLOAT, "stat_sou", "Source statics [ms]" );
static csHeaderInfo const HDR_STAT_TOT     ( TYPE_FLOAT, "stat_tot", "Total statics [ms]" );
static csHeaderInfo const HDR_STAT_RES     ( TYPE_FLOAT, "stat_res", "Residual statics [ms]" );

static csHeaderInfo const HDR_TIME_CODE     ( TYPE_INT, "time_code",  "Time code (1: local, 2: GMT, 3: Other, 4: UTC)" );
static csHeaderInfo const HDR_TIME_YEAR     ( TYPE_INT, "time_year",  "Year of shot/record" );
static csHeaderInfo const HDR_TIME_DAY      ( TYPE_INT, "time_day",   "Julian day of shot/record" );
static csHeaderInfo const HDR_TIME_HOUR     ( TYPE_INT, "time_hour",  "Hour of shot/record" );
static csHeaderInfo const HDR_TIME_MIN      ( TYPE_INT, "time_min",   "Minute of shot/record" );
static csHeaderInfo const HDR_TIME_SEC      ( TYPE_INT, "time_sec",   "Second of shot/record" );
static csHeaderInfo const HDR_TIME_MSEC     ( TYPE_INT, "time_msec",  "Millisecond shot/record time fraction [ms]" );
static csHeaderInfo const HDR_TIME_USEC     ( TYPE_INT, "time_usec",  "Microsecond shot/record time fraction [us]" );
static csHeaderInfo const HDR_TIME_NANO     ( TYPE_INT, "time_nano",  "Nanosecond shot/record time fraction [ns]" );
static csHeaderInfo const HDR_TIME_SAMP1    ( TYPE_INT, "time_samp1",    "Absolute time of first sample [s], seconds since 01-Jan-1970" );
static csHeaderInfo const HDR_TIME_SAMP1_US ( TYPE_INT, "time_samp1_us", "Absolute time of first sample [us], fraction" );
static csHeaderInfo const HDR_TAPENO        ( TYPE_INT, "tapeno", "Tape number" );
static csHeaderInfo const HDR_TRC_TYPE      ( TYPE_INT, "trc_type", "Trace type" );
static csHeaderInfo const HDR_TRCNO         ( TYPE_INT, "trcno", "Trace number" );
//static csHeaderInfo const HDR_TRCSEQ        ( TYPE_INT, "trcseq", "Trace sequence number" );
//static csHeaderInfo const HDR_TRCSEQ_LINE   ( TYPE_INT, "trcseq_line", "Trace sequence number within line" );
static csHeaderInfo const HDR_TRC_EDIT      ( TYPE_INT, "trc_edit", "Trace edit code" );

static csHeaderInfo const HDR_UNIT_COORD    ( TYPE_INT, "unit_coord", "Coordinate unit (1: length, 2: arc sec, 3: decimal deg, 4: deg, min, sec)" );

static csHeaderInfo const HDR_WDEP          ( TYPE_FLOAT, "wdep", "Water depth [m]" );

template <typename T> class csVector;

class csStandardHeaders {
public:
  /**
  * @param name (i) Name of standard header
  * @return true if standard header with given name exists
  */
  static bool isStandardHeader( std::string const& name );
  /**
  * @param info (i) Header info of standard header
  * @return true if standard header with given exists that matches given info name, type and description
  */
  static bool isStandardHeader( csHeaderInfo const* info_in );
  /**
  * Get standard header info
  * @param name (i) name of standard header
  * @return Constant pointer to header info of standard header. Returns NULL if no standard header with given name was found
  */
  static csHeaderInfo const* get( std::string const& name );
  /**
  * Get all standard headers
  * @param hdrList (o) List of header infos of all existing standard headers
  */
  static void getAll( csVector<csHeaderInfo const*>* hdrList );

 private:
//-------------------------------------------------
//  template <typename T> static int findPosition( T const* value, T const* values, int nValues ) {
  /**
  * Helper method
  * Find name in given header info list and return position index
  * @param name     (i) Name of header to search
  * @param infoList (i) List of header infos that shall be searched. This list must be sorted by name!
  * @return position index of given name. Returns -1 if given name is 'smaller' than first header name, nValues if name is 'larger' than all header names.
  */
  static int findPosition( std::string const& name, csVector<csHeaderInfo const*> const* infoList );

};

} // namespace

#endif

